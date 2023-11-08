from sklearn.feature_selection._mutual_info import _compute_mi_cc
from itertools import combinations
import numpy as np

def get_cluster_mean_mi(members,data,n_neighbors=8):
    mis = []
    for x,y in combinations(members,2):
        x_data = data.loc[x].values
        y_data = data.loc[y].values
        mis.append(_compute_mi_cc(x_data,y_data,n_neighbors))
    return np.mean(mis)
