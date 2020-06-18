# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/20
# content:    Test the algorithm on same artificial data
import numpy as np
import pandas as pd
import anndata
import anndata_utils


def test_split():
    X = np.array([
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 4],
        [3, 5],
        ])
    obs = pd.DataFrame([
        ['blue'],
        ['yellow'],
        ['blue'],
        ['blue'],
        ['blue'],
        ],
        columns=['color'])
    adata = anndata.AnnData(X=X, obs=obs)
    adata.var_names = ['Gene1', 'Gene2']

    dic = anndata_utils.split(adata, 'color', axis='obs')

    assert(list(dic.keys()) == ['blue', 'yellow'])
