"""
contains code that combines all steps for a complete analysis
"""
from . import mixture_modelling as mm
from . import bootstrap as bs
from . import process_normalise as pn
from scipy.stats import zscore
import multiprocessing as mp


def main(data,clust_ratio,annot_df=None,
    covar_type='diag',n_init=5,fit_fn=None,n_bootstraps=500,
    bs_subsample_size=None,bs_subsample_replacement=False,
    bs_processes=1,bs_score_fn=None,bs_mem_fn=None,
    min_size=2,min_maxloc=9,mi_thresh=0,seed=None,
    clusttable_fn = None, membertable_fn = None, pdf_fn=None):

    """
    Main function to run full GIP clustering workflow on dataset 

    Args:
        data (pd.DataFrame): 
            input data to perform resampling and clustering on
            "raw" protein abundances intensities per fraction
            index should contain protein identifiers
        clust_ratio (float): to determine the number of output clusters
            ratio of the number of clusters to proteins
            a ratio of 0.5 would mean one cluster per two detected 
            proteins in the dataset
        annot_df (pd.DataFrame, optional): dataframe with protein annotations
            row index values should match protein ids in "data"
            columns contain additional annotations of proteins
            that will be included in the GIP output
        covar_type (string, optional): defaults to 'diag'
            the covar_type used in GMM clustering, as described
            in sklearn.mixture.GaussianMixture
        n_init (int, optional): defaults to 5
            the number of initializations used in GMM clustering,
            as described in sklearn.mixture.GaussianMixture
        fit_fn (str, optional): 
            filepath to store pickle of initial clustering fit
            if the file already exists it will load the fit from
            this file rather than running initial clustering
        n_bootstraps (int, optional):
            number of bootstrap iterations to perform
        bs_subsample_size (int, optional): Defaults to None.
            Size of subsample taken when resampling
        bs_subsample_replacement (bool, optional): Default False (without).
            to perform resampling with or without replacement 
        bs_processes (int, optional):Defaults to 1.
            number of threads/processes for bootstrapping procedure
        bs_score_fn (string, optional): Defaults to None.
            filename to store or retrieve bootstrap scores
            if both score_fn and membership_fn files already exist:
                results retrieved from file, bootstrapping skipped
            if file does not exist yet:
                bootstrapping results are stored in these files
        bs_mem_fn (string, optional): Defaults to None.
            filename to store or retrieve bootstrapped cluster membership
            if both score_fn and membership_fn files already exist:
                results retrieved from file, bootstrapping skipped
            if file does not exist yet:
                bootstrapping results are stored in these files
        min_size (int, optional): Defaults to 2
            minimum cluster size. any smaller clusters are filtered out
        min_maxloc (int, optional): defaults to 9
            clusters with its max peak at fraction number lower than
            'min_maxloc' are filtered out
        mi_thresh (float, optional): defaults to 0
            clusters with a average mutual information lower than
            'mi_thresh' are filtered out
        seed (int, optional): Defaults to None.
            random seed used for resampling and clustering
        clusttable_fn (str, optional):
            if provided will save table with output clusters to given
            filepath
        membertable_fn (str, optional):
            if provided will save table with output cluster members to
            given filepath
        pdf_fn (str, optional):
            if provided will save a pdf with plots for all output 
            clusters to given filepath.

    Returns: dict with following entries
        clust_res (dict): dict with initial clustering results,
        bs_res (dict): dict with bootstrapping results,
        clusttable (pd.DataFrame): 
            output table with overview of result clusters
        membertable (pd.DataFrame): 
            output table with overview of clustered proteins
        cluster_flows (dict):
            contains degree of member overlap between clusters during
            bootstrapping procedure, for all clusters that have 
            overlapping cluster members during bootstrapping
    """

    # set multiprocessing method so mp works with sklearn models
    try:
        mp.set_start_method('forkserver')
    except RuntimeError:
        print(('failed to set forkserver multiprocessing method, '
             'assuming it was already set'))
    except Exception as e:
        msg = f'unexpected error setting multiprocessing method: {e}'
        raise ValueError(msg)

    # scale data, get protein abundances
    scaled_data = pn.z_score(data)
    relative_abuns = pn.get_relative_abuns(data)

    # perform GMM cluster analysis
    print('performing original clustering run..')
    clust_res = mm.mixture_modelling(scaled_data,clust_ratio=clust_ratio,
                        covar_type=covar_type,
                        n_init=n_init,fit_fn=fit_fn,seed=seed)
    
    # perform bootstrapping
    print('performing bootstrap analysis..')
    bs_res = bs.run_bootstrap(scaled_data,clust_res['pred_dict'],
        n=n_bootstraps,processes=bs_processes,
        subsample_size=bs_subsample_size,
        replacement=bs_subsample_replacement,
        score_fn=bs_score_fn,membership_fn=bs_mem_fn,seed=seed)

    # filter clusters
    print('processing cluster and bootstrapping results..')
    clusttable = mm.filter_clusters(
        clust_res['feature_df'],min_size=min_size,
        min_maxloc=min_maxloc)

    # compute mutual information
    mean_mis = mm.compute_cluster_mean_mis(
        clusttable.index.values,
        clust_res['pred_dict'],scaled_data)

    # compute mean relative abundances
    mean_abuns = pn.cluster_mean_abuns(relative_abuns,clust_res['pred_dict'])

    # cluster table
    clusttable['mean_mi'] = clusttable.index.map(mean_mis)
    clusttable['stability'] = clusttable.index.map(bs_res['stabilities'])
    clusttable['mean_abundance'] = clusttable.index.map(mean_abuns)
    clusttable['mean_abundance_z'] = zscore(clusttable['mean_abundance'])
    clusttable.index.name='clust_id'

    # filter on mutual information of 0
    if mi_thresh != None:
        print(('filtering out clusters with mutual'
         f' information <= {mi_thresh}'))
        clusttable = clusttable[clusttable['mean_mi']>mi_thresh]

    # clustmember table
    membertable = mm.create_clustmember_table(
        clust_res['pred'],annot_df=annot_df)
    membertable = mm.filter_clustmem_table(membertable,
                                    clusttable.index.values)
    membertable['rel_abundance'] = membertable['identifier'].map(
                                                        relative_abuns.to_dict())
    membertable['rel_abundance_z'] = zscore(membertable['rel_abundance'])
    membertable['frequency'] = membertable['identifier'].map(bs_res['flat_freqs'])
    # set any proteins missing without freq value to 0
    membertable.loc[:,'frequency'].fillna(0,inplace=True)
    
    # compute cluster flows, only for clusters that passed filtering
    # sort filtered ids by clust size for their use in generating pdf of plots
    filtered_ids = clusttable.sort_values(by='size',ascending=False).index.values
    filt_clusters = {key:val for key,val in clust_res['pred_dict'].items()
                if key in filtered_ids}
    cluster_flows = bs.compute_cluster_flows(filt_clusters,bs_res['memberships'])

    # write the cluster and membertables to file
    if membertable_fn:
        membertable.to_csv(membertable_fn, sep='\t')
    if clusttable_fn:
        clusttable.to_csv(clusttable_fn, sep='\t')

    # generate pdf with plots of all clusters that passed filtering
    if pdf_fn:
        print(f'generating pdf with plotted clusters: {pdf_fn}')
        means = clust_res['fit'].means_
        covs = clust_res['fit'].covariances_
        mm.generate_plot_pdf(pdf_fn,filtered_ids,means,covs,
            clust_res['pred'],scaled_data)

    return {
        'clust_res':clust_res,
        'bs_res':bs_res,
        'clusttable':clusttable,
        'membertable':membertable,
        'cluster_flows':cluster_flows,
    }
