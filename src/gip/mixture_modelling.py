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
    get locations of maximum of each cluster
    """
    return np.argmax(means,axis=1)

def get_sizes(labels):
    """
    get cluster sizes from predicted assignments
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
    return pd.DataFrame([sizes,maxlocs],index = ['size','maxloc']).T

def filter_clusters(feature_df,min_size,min_maxloc):
    size_filt = feature_df[feature_df['size']>=min_size]
    peak_filt = size_filt[size_filt['maxloc']>=min_maxloc]
    return peak_filt.copy()

def save_fit(fit,filename):
    pickle.dump(fit,open(filename,'wb'))

def load_fit(filename):
    return pickle.load(open(filename,'rb'))

def mixture_modelling(data,clust_ratio,covar_type='diag',
                        n_init=5,fit_fn=None,seed=None):
    """
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

def mixture_modelling_with_filt(data,clust_ratio,covar_type='diag',
    n_init=5,min_size=2,min_maxloc=9,
    fit_fn = None,seed=None):
    """
    """
    # fit and predict
    if fit_fn:
        fit = load_fit(fit_fn)
    else:
        n_clusts = round(data.shape[0]*clust_ratio)
        fit = GaussianMixture(n_clusts,n_init=n_init,
                covariance_type=covar_type,random_state=seed).fit(data)
    
    pred = pd.Series(fit.predict(data),index=data.index)
    pred_dict = labels_to_dict(pred)

    # get clust properties and filter
    maxlocs = find_maxlocs(fit.means_)
    sizes = get_sizes(pred)
    feature_df = get_feature_df(sizes,maxlocs)
    filtered_clusters = filter_clusters(feature_df,min_size,min_maxloc)
    
    mean_mis = compute_cluster_mean_mis(
        filtered_clusters.index.values,pred_dict,data)
    
    filtered_clusters['mean_mi'] = filtered_clusters.index.map(mean_mis)
    
    # return results
    return {
        'fit':fit,
        'pred':pred,
        'pred_dict':pred_dict,
        'feature_df':feature_df,
        'filtered_clusters':filtered_clusters,
    }

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