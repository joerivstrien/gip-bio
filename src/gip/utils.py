def labels_to_dict(labels):
    """
    loads ComplexSet from labels

    labels: pd.Series. index are proteins,
            values are cluster labels
    """
    clusters = labels.unique()
    as_dict = {clust:labels[labels==clust].index.to_list()
               for clust in clusters}
    return as_dict

def write_gmt(complexes,fn):
    """"""
    with open(fn,'w') as f_obj:
        for name,members in complexes.items():
            memberstring = "\t".join(members)
            f_obj.write(f'{name}\tnan\t{memberstring}\n')
