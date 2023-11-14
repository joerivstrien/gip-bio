# base library imports
import pickle
import matplotlib.backends.backend_pdf
from matplotlib import pyplot as plt
from os.path import exists

# third pary imports

import pandas as pd
import numpy as np
import seaborn as sns
sns.set_style('whitegrid')

from sklearn.mixture import GaussianMixture

# local imports
from .mutual_information import get_cluster_mean_mi
from .utils import labels_to_dict

def find_maxlocs(means):
    """
    get locations of maximum of each cluster's mean migration pattern

    Args:
        means (2d numpy array): mean migration pattern of each cluster

    Returns:
        1d numpy array: location of migration pattern maximum 
                        for each cluster
    """
    return np.argmax(means,axis=1)

def get_sizes(labels):
    """
    get cluster sizes from predicted assignments

    Args:
        labels (pd.Series): cluster label for each clustered protein

    Returns:
        numpy array: size of each cluster
    """
    return labels.value_counts(sort=False).sort_index().values

def draw_cluster(ax,mean,std):
    """
    draws mean +- std of cluster
    """
    x = list(range(len(mean)))
    ax.plot(mean,linewidth=1)
    ax.fill_between(x,mean-2*std, mean+2*std,alpha=0.7)

def generate_plot_pdf(outfn,clust_ids,means,stds,labels,abun_data):
    """
    Generates pdf with migration plot for each cluster

    Args:
        outfn (string): filepath of output pdf file 
        clust_ids (list): list of cluster ids to plot
        means (numpy array): 
            mean values of each cluster for each fraction
        stds (numpy array): 
            covariance values for each cluster for each fraction
        labels (dict): cluster assignments for each protein
        abun_data (pd.DataFrame):
            scaled complexome profiling data used as input for clustering
    """
    pdf = matplotlib.backends.backend_pdf.PdfPages(outfn)

    for clust in clust_ids:
        fig,ax = plt.subplots()
        draw_cluster(ax,means[clust],stds[clust])

        members = np.where(labels==clust)[0]

        prot_data = abun_data.values[members,:]
        for prot in prot_data:
            ax.plot(prot, color='grey',alpha=0.5)
        ax.set_title(f'cluster: {clust}, size: {len(members)}')

        pdf.savefig(fig)

    pdf.close()

def get_feature_df(sizes,maxlocs):
    """
    generate dataframe with cluster feature annotations

    Args:
        sizes (numpy array): size of each cluster
        maxlocs (numpy array): 
            location of migration pattern maximum 
            for each cluster

    Returns:
        pd.DataFrame: table with row for each cluster
            columns contain size and maxloc for each cluster 
    """
    return pd.DataFrame([sizes,maxlocs],index = ['size','maxloc']).T

def mixture_modelling(data,clust_ratio,covar_type='diag',
                        n_init=5,fit_fn=None,seed=None):
    """
    Perform Gaussian mixture modeling clustering on complexome profile

    Args:
        data (pd.DataFrame): 
            input data to perform resampling and clustering on
            "raw" protein abundances intensities per fraction
            index should contain protein identifiers
        clust_ratio (float): to determine the number of output clusters
            ratio of the number of clusters to proteins
            a ratio of 0.5 would mean one cluster per two detected 
            proteins in the dataset
        covar_type (string, optional):
            the covar_type used in GMM clustering, as described
            in sklearn.mixture.GaussianMixture
        n_init (int, optional): defaults to 5
            the number of initializations used in GMM clustering,
            as described in sklearn.mixture.GaussianMixture
        fit_fn (str, optional): 
            filepath to store pickle of initial clustering fit
            if the file already exists it will load the fit from
            this file rather than running initial clustering
        seed (int, optional): Defaults to None.
            random seed used for resampling and clustering

    Returns: dict with following entries:
        fit (sklearn.mixture.GaussianMixture): the fitted mixture model 
        pred (pd.Series): the predicted cluster assignment of each protein 
        pred_dict (dict): the predicted cluysters as a dictionary
            keys: cluster ids, items: lists with protein members 
        feature_df (pd.DataFrame): table with cluster sizes and maxlocs 
    """

    # if no fit_fn is provided: run model and dont save fit
    if not fit_fn:
        n_clusts = round(data.shape[0]*clust_ratio)
        fit = GaussianMixture(
            n_clusts,n_init=n_init,random_state=seed,
            covariance_type=covar_type).fit(data)
    else:
        # if fit_fn and file already exists: load fit from file
        if exists(fit_fn):
            fit = load_fit(fit_fn)
        # if fit_fn but file does not exist: run model and save fit
        else:
            n_clusts = round(data.shape[0]*clust_ratio)
            fit = GaussianMixture(
                n_clusts,n_init=n_init,random_state=seed,
                covariance_type=covar_type).fit(data)
            save_fit(fit,fit_fn)

    pred = pd.Series(fit.predict(data),index=data.index)
    pred_dict = labels_to_dict(pred)

    # get clust properties
    maxlocs = find_maxlocs(fit.means_)
    sizes = get_sizes(pred)
    feature_df = get_feature_df(sizes,maxlocs)

    # return results
    return {
        'fit':fit,
        'pred':pred,
        'pred_dict':pred_dict,
        'feature_df':feature_df,        
    }

def filter_clusters(feature_df,min_size,min_maxloc):
    """
    filter cluster dataframe based on cluster size and maxloc

    Args:
        feature_df (pd.DataFrame): clusters (rows) with features (columns)
        min_size (int): clusters with lower size are removed
        min_maxloc (int): clusters with lower maxloc are removed

    Returns:
        pd.DataFrame: feature_df with filtered out clusters removed
    """
    size_filt = feature_df[feature_df['size']>=min_size]
    peak_filt = size_filt[size_filt['maxloc']>=min_maxloc]
    return peak_filt.copy()

def save_fit(fit,filename):
    pickle.dump(fit,open(filename,'wb'))

def load_fit(filename):
    return pickle.load(open(filename,'rb'))

def create_clustmember_table(labels,annot_df=None):
    membertable = pd.DataFrame(labels)
    membertable.reset_index(inplace=True)
    membertable.columns = ['identifier','clust_id']
    membertable.set_index('clust_id',inplace=True)

    if isinstance(annot_df,pd.DataFrame):
        membertable = membertable.merge(annot_df,left_on='identifier',
            right_index=True,how='left')

    return membertable

def filter_clustmem_table(clustmem_table,filtered_ids):
    return clustmem_table.loc[filtered_ids]

def compute_cluster_mean_mis(clust_ids,member_dict,data,neighbors=8):
    clust_mean_mis = {}
    for cid in clust_ids:
        members = member_dict[cid]
        mean_mi = get_cluster_mean_mi(members,data,neighbors)
        clust_mean_mis[cid] = mean_mi
    return clust_mean_mis