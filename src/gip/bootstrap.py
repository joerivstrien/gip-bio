import pandas as pd
from numpy.random import default_rng
from sklearn.mixture import GaussianMixture
from .utils import labels_to_dict
import multiprocessing as mp
import json
from os.path import exists

"""
BOOTSTRAPPING PROCEDURE
"""

def run_bootstrap(data,original_clusters,subsample_size=None,n=500,
                  replacement=False,
                  processes=1,score_fn=None,membership_fn=None,
                  seed=None):
    """
    perform bootstrapping workflow, repeated reclustering of data

    Args:
        data (pd.DataFrame): input data to perform resampling and clustering on
        original_clusters (dict of lists): 
            clusters resulting from original clustering
            each entry (cluster) contains list with protein ids (string)
        subsample_size (int, optional): Defaults to None.
            Size of subsample taken when resampling
        n (int, optional): Defaults to 500.
            number of bootstrap iterations
        replacement (bool, optional): Defaults to False (without).
            to perform resampling with or without replacement
        processes (int, optional):Defaults to 1.
            number of threads/processes
        score_fn (string, optional): Defaults to None.
            filename to store or retrieve bootstrap scores
            if both score_fn and membership_fn files already exist:
                results retrieved from file, bootstrapping skipped
            if file does not exist yet:
                bootstrapping results are stored in these files
        membership_fn (string, optional): Defaults to None.
            filename to store or retrieve bootstrapped cluster membership
            if both score_fn and membership_fn files already exist:
                results retrieved from file, bootstrapping skipped
            if file does not exist yet:
                bootstrapping results are stored in these files
        seed (int, optional): Defaults to None.
            random seed used for resampling and clustering

    Returns:
        dict: bootstrapping results
            scores: pd.DataFrame
                overlap scores of each original cluster with best
                fitting bootstrapped cluster, for each iteration
                columns are orig clusters, rows are iterations
            memberships: dict of dicts
                frequency with which each protein is part of each cluster
                each outer entry is a cluster
                nested dicts are each member with associated frequency                
            stabilities: dict
                stores cluster stability (float, mean overlap) per cluster
            flat_freqs: dict
                stores the frequency of each protein in their original cluster
    """
    if (not score_fn) or (not membership_fn):
        print('performing bootstrap')
        scores,memberships = bootstrap_cluster_result(
            data,original_clusters,
            subsample_size=subsample_size,n=n,
            replacement=replacement,
            processes=processes,seed=seed)

        if score_fn:
            save_bootstrap_scores(scores,score_fn)
        if membership_fn:
            save_membership(memberships,membership_fn)

    else:
        if exists(score_fn) and exists(membership_fn):
            print('bootstrap results available. loading from file..')
            scores = load_bootstrap_scores(score_fn)
            memberships = load_membership(membership_fn)

        else:
            print('performing bootstrap')
            scores,memberships = bootstrap_cluster_result(
                data,original_clusters,
                subsample_size=subsample_size,n=n,
                replacement=replacement,
                processes=processes,seed=seed)
            save_bootstrap_scores(scores,score_fn)
            save_membership(memberships,membership_fn)

    stabilities,flat_freqs = process_bootstrap_results(
        scores,memberships,original_clusters
    )

    return {
        'scores':scores,
        'memberships':memberships,
        'stabilities':stabilities,
        'flat_freqs':flat_freqs,
    }

def process_bootstrap_results(scores,memberships,original_clusters):
    """
    processes bootstrap results to metrics for each cluster and protein


    Returns:
        stabilities: dict
            stores cluster stability (float, mean overlap) per cluster
        flat_freqs: dict
            stores the frequency of each protein in their original cluster
    """

    # compute cluster stabilities from bootstrap scores
    stabilities = scores.mean().to_dict()
        
    # get flat frequencies for only original cluster members
    flat_frequencies = {}
    for cid,member_freqs in memberships.items():
        for mem,freq in member_freqs.items():
            if mem in original_clusters[cid]:
                flat_frequencies[mem] = freq

    return stabilities,flat_frequencies
  
def bootstrap_cluster_result(data,original_clusters, subsample_size=None,
                             replacement=False,n=500,processes=1,seed=None):
    """
    perform bootstrapped reclustering on initial clustering and data

    see run_bootstrap documentation for parameter descriptions

    Returns:
        tuple: (scores, memberships)
            scores: pd.DataFrame
                overlap scores of each original cluster with best
                fitting bootstrapped cluster, for each iteration
                columns are orig clusters, rows are iterations
            memberships: dict of dicts
                frequency with which each protein is part of each cluster
                each outer entry is a cluster
                nested dicts are each member with associated frequency                
    """
    rng = default_rng(seed=seed)
    n_clusts = len(original_clusters)

    if not subsample_size:
        subsample_size = round(data.shape[1]/2)

    to_iter = ((data,original_clusters,n_clusts,subsample_size,
                replacement,rng,seed) for _ in range(n))
    pool = mp.Pool(processes)
    scores_members = pool.starmap(perform_bootstrap, to_iter)

    #unzip scores_members
    scores, members = zip(*scores_members)

    # process scores
    score_df = pd.DataFrame(list(scores))

    # process member counts
    membership_counts = get_membership_counts(list(members),original_clusters)
    normed_memcounts = normalize_memcounts(membership_counts,n)

    return score_df, normed_memcounts

def perform_bootstrap(data,original_clusters,n_clusts,
                      subsample_size,replacement,rng,seed):

    print('starting (another) bootstrap round..')
    # resample dataset
    subsample = resample_data(data,subsample_size,replacement,rng)
    
    # perform clustering
    bootstrap_clusters = run_GMM(subsample,n_clusts,seed=seed)

    # determine max jaccards for original clusters
    bootstrap_scores,boostrap_members = jaccard_target_query(
        original_clusters,bootstrap_clusters)
    
    return bootstrap_scores,boostrap_members

def run_GMM(data,n_clusts,seed=None):
    """
    function to run GMM, returns dict of clusters.

    'diag' covar type is hardcoded
    """
    fit = GaussianMixture(n_clusts,covariance_type='diag',
                          random_state=seed).fit(data)
    pred = pd.Series(fit.predict(data),index=data.index)
    pred_dict = labels_to_dict(pred)
    return pred_dict    

def save_bootstrap_scores(bootstrap_scores,fname):
    bootstrap_scores.to_csv(fname,sep='\t',index=False)

def load_bootstrap_scores(fname):
    df = pd.read_csv(fname,sep='\t')
    df.columns = df.columns.astype(int)
    return df

def resample_data(data, subsample_size,replacement,rng):
    """
    taking a subsample of the data (columns), of given size
    """
    choice = rng.choice(data.columns.values,subsample_size,
                        replace=replacement)
    return data.loc[:,choice]

"""
CODE RELATED TO BOOTSTRAP MEMBERSHIP COUNTING
"""

def get_membership_counts(bootstrap_members,original_clusters):
    """
    """
    membership_counts = {cid:dict() for cid in original_clusters.keys()}

    # loop over bootstraps
    for bclust_members in bootstrap_members:
        for clust in original_clusters.keys():
            for member in bclust_members[clust]:
                if member in membership_counts[clust]:
                    membership_counts[clust][member] += 1
                else:
                    membership_counts[clust][member] = 1

    return membership_counts

def normalize_memcounts(membership_counts,n_bootstraps):
    """
    """
    normed_counts = {}
    for cid, clust_dict in membership_counts.items():
        norm_clust_dict = {}
        for mid,count in clust_dict.items():
            norm_clust_dict[mid] = clust_dict[mid]/n_bootstraps
        normed_counts[cid] = norm_clust_dict
    return normed_counts

def save_membership(membership_counts,fname):
    membership_counts = {int(key):val for key,val in membership_counts.items()}
    with open(fname,'w') as f_obj:
        f_obj.write(json.dumps(membership_counts))    

def load_membership(fname):
    with open(fname) as json_file:
        memberships = json.load(json_file)
        return {int(key):val for key,val in memberships.items()}

"""
code to compute bootstrapped "cluster flow"
"""
def compute_cluster_flows(clusters,member_weights):
    """
    computes cluster flow between overlapping clusters
    """
    # determine which clusters have any overlap
    overlapping = get_overlapping_clusters(clusters,member_weights)

    # computing flow for overlapping clusters
    cluster_flows = {}
    for left,right in overlapping:
        flow = compute_flow(clusters[left],member_weights[right],
                            clusters[right],member_weights[right])
        cluster_flows[frozenset([left,right])] = flow
    cluster_flows = pd.Series(cluster_flows)

    return cluster_flows


def get_overlapping_clusters(original_clusters,normalised_counts):
    """
    """
    overlapping_clusters = set([])

    for cid,members in original_clusters.items():
        to_compare = set(original_clusters.keys())
        to_compare.remove(cid)
        for comp_id in to_compare:
            comp_members = normalised_counts[comp_id]
            if set(comp_members.keys()) & set(members):
                overlapping_clusters.update({frozenset([cid,comp_id])})

    return overlapping_clusters

def compute_flow(left_real,left_weights,right_real,right_weights):
    """
    compute the "bootstrap flow" between two bootstrapped clusters
    [left/right]_real: list of cluster members
    [left/right]_weights: dict with members and their bootstrapped weights
    """

    # get cluster flow
    left_to_right = sum([right_weights[mid] for mid in left_real if mid in right_weights])
    right_to_left = sum([left_weights[mid] for mid in right_real if mid in left_weights])

    # get cluster weight (should be only of real members)
    # left_weight = sum([left_weights[mid] if mid in left_weights else 0 for mid in left_real])
    # right_weight = sum([right_weights[mid] if mid in right_weights else 0 for mid in right_real])

    # compute flow "quotient?" think about a name
    flow = (left_to_right + right_to_left)/(len(left_real) + len(right_real))
    return flow

"""
JACCARD COEFFICIENT COMPUTATION
"""

def jaccard_target_query(query,target):
    """
    computes max jaccard coef for each query cluster against all target clusters
    """
    max_scores = {}
    max_members = {}
    for name,members in query.items():
        max_score,members = jaccard_one_many(members,target)
        max_scores[name] = max_score
        max_members[name] = members
    return max_scores,max_members

def jaccard_one_many(one,many):
    """
    returns best match jaccard score and cluster members

    DOES NOT HANDLE CASE WHEN MORE THAN ONE ARE TIED FOR MAX YET
    """
    scores = []
    cur_max = 0
    cur_max_members = []
    for members in many.values():
        score = jaccard_coefficient(one,members)
        if score > cur_max:
            cur_max = score
            cur_max_members = members
    
    return cur_max,cur_max_members

def jaccard_coefficient(A,B):
    """"""
    intersect = set(A) & set(B)
    union = set(A) | set(B)
    if len(union) == 0:
        return 0
    return len(intersect)/len(union)



if __name__ == "__main__":
    original_clusters = {
        1:['a','b'],
        2:['d','e'],
    }
    bootstrap_members = [
        {1:['a'],2:['d','e']},
        {1:['a','b','d'],2:['e']},
    ]

    mem_counts = get_membership_counts(bootstrap_members,original_clusters)
    normed_counts = normalize_memcounts(mem_counts,2)
    print('normalised bootstrap counts')
    print(normed_counts)

    print('original clusters')
    print(original_clusters)

    get_overlapping_clusters(original_clusters,normed_counts)

    # compute flow between two bootstrapped clusters
    flow = compute_flow(original_clusters[1],normed_counts[1],
                 original_clusters[2],normed_counts[2])

    print(flow)

    quit()

    # A = [1,2,3]
    # B = [3,4,5]
    # print(jaccard_coefficient(A,B))

    # one = [1,2,3]

    # many = {
    #     1: {1,2},
    #     2: {3},
    #     3: {4,5,6}
    #     }

    # print(jaccard_one_many(one,many))

    # target = {
    #     1: {1,2,6},
    #     2: {7,8},
    #     3: {9,3}
    #     }

    # query = {
    #     1: {1,2,3},
    #     2: {7,8,9},
    #     3: {6}
    #     }

    # print(jaccard_target_query(query,target))
    from mixture_modelling import *
    mp.set_start_method('forkserver')

    """TRY MULTIPROCESSING"""
    annot_fn = '../data/bigsearch_human_plasmo_annot.tsv'
    annot_df = pd.read_csv(annot_fn,sep='\t',index_col=0)

    abs_all_fn = '../data/profiles/pf_m_abs200_HR_abuns.tsv'
    abs_all_df = pd.read_csv(abs_all_fn,sep='\t',index_col=0)
    abs_all_z = z_score(abs_all_df)

    abs_all_results = mixture_modelling(abs_all_z,0.6,
                                    fit_fn='../results/abs200_all_fit_06_13oct.pickle')
    real_clusters = abs_all_results['pred_dict']

    n_clusters = round(abs_all_z.shape[0]*0.6)
    bootstrap_res = bootstrap_cluster_result_mp(abs_all_z,real_clusters,n_clusters,processes=6,
                                                n=10)