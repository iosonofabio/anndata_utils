# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/20
# content:    Kolmogorov Smirnov test on gene expression for AnnData objects
import numpy as np
import pandas as pd


def split(adata, columns, axis='obs'):
    '''Split AnnData by metadata column

    Args:
        adata (AnnData): The object to split
        columns (str or list): column(s) of metadata table to split by
        axis (str): 'obs' (default) or 'var'

    Returns:
        dict of AnnData with keys equal to the unique values found in that
        columns of the appropriate metadata (adata.obs or adata.var)

    '''
    if axis not in ('obs', 'var'):
        raise ValueError('axis must be "obs" or "var"')

    meta = getattr(adata, axis)

    if isinstance(columns, str):
        column = columns
        if column not in meta.columns:
            raise ValueError('column not found')
        metac = meta[column]
    else:
        if len(columns) < 1:
            return adata
        elif len(columns) == 1:
            return split(adata, columns[0], axis=axis)

        # Merge with @ (a hack, but alright in most cases)
        metac = meta[columns[0]].astype(str)
        for i, col in enumerate(columns[1:]):
            metac = metac + '@' + meta[col].astype(str)

    # Unique values
    metau = metac.unique()

    # Construct dict
    d = {}
    for key in metau:
        idx = metac.index[metac == key]
        if axis == 'obs':
            asub = adata[idx]
        else:
            asub = adata[:, idx]

        if not isinstance(columns, str) and (len(columns) > 1):
            key = key.split('@')
            for i, col in enumerate(columns):
                key[i] = type(meta[col].iloc[0])(key[i])
            key = tuple(key)

        d[key] = asub

    return d


def expressing_fractions(adata, columns, axis='obs', greater_than=0):
    '''Fraction of expressors by metadata column

    Args:
        adata (AnnData): The object to split
        columns (str or list): column(s) of metadata table to split by
        axis (str): 'obs' (default) or 'var'
        greater_than (float): only expressors stricly above this threshold are
          counted

    Returns:
        pandas DataFrame with rows equal to the non-integrated axis names
        and columns equal to the keys of the AnnData split.

    '''
    adatad = split(adata, columns, axis=axis)

    if axis == 'var':
        iax = 1
        index = adata.obs_names
    else:
        iax = 0
        index = adata.var_names

    fracd = {}
    for key, adatai in adatad.items():
        fr = np.asarray((adatai.X > greater_than).mean(axis=iax))[0]
        fracd[key] = fr
    fracd = pd.DataFrame(fracd, index=index)

    return fracd

