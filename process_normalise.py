"""
for now this is where I do data processing, normalisation etc.
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler


# load data
def parse_profile(tsv_fn):
    """
    parses standard format complexome profile into df

    Args:
        tsv_fn (str): filepath of file containing:
            complexome profile in tab separated text format.
            single header col with fraction ids
            single index row with protein ids
            numeric abundance values
    Returns:
        pd df: complexome profile as dataframe
    """
    df = pd.read_csv(tsv_fn, sep='\t')
    df.iloc[:, 0] = df.iloc[:, 0].astype(str)
    df.set_index(df.columns[0], inplace=True)
    df.index.name = 'prot_ids'
    if df.empty:
        msg = f'parsed profile contains no values: {tsv_fn}'
        raise ValueError(msg)
    try:
        df = df.astype(float)
    except BaseException:
        msg = f'profile contains non-numeric values: {tsv_fn}'
        raise ValueError(msg)
    return df

# z-score normalisation
# from Iris 
def z_score(df):
    """
    perform z-score normalization: (x-avg(x))/std -> mean = 0, std = 1

    Args:
        df (dataframe): input data matrix with slices as columns and protein ids as rows
        row_names (list): protein ids
    """
    z = StandardScaler().fit_transform(df.T)
    z = pd.DataFrame(z.T,index=df.index)
    return z

# get relative protein abundances
def get_relative_abuns(profile,log=True):
    abuns = profile.sum(axis=1)
    relative_abuns = abuns/abuns.sum()

    if log:
        return np.log(relative_abuns)
    else:
        return relative_abuns

def cluster_mean_abuns(protein_abuns,clusters):
    mean_abuns = {}
    for cid,members in clusters.items():
        mean = protein_abuns.loc[members].mean()
        mean_abuns[cid] = mean
    return mean_abuns


