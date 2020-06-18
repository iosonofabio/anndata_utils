# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/20
# content:    Kolmogorov Smirnov test on gene expression for AnnData objects
import numpy as np
import pandas as pd


def split(adata, column, axis='obs'):
    '''Split AnnData by metadata column

    Args:
        adata (AnnData): The object to split
        column (str): column of metadata table to split by
        axis (str): 'obs' (default) or 'var'

    Returns:
        dict of AnnData with keys equal to the unique values found in that
        columns of the appropriate metadata (adata.obs or adata.var)

    '''
    if axis not in ('obs', 'var'):
        raise ValueError('axis must be "obs" or "var"')

    meta = getattr(adata, axis)
    if column not in meta.columns:
        raise ValueError('column not found')
    meta = meta[column]

    # Unique values
    metau = meta.unique()

    # Construct dict
    d = {}
    for key in metau:
        idx = meta.index[meta == key]
        if axis == 'obs':
            asub = adata[idx]
        else:
            asub = adata[:, idx]

        d[key] = asub

    return d
