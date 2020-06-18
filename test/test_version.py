# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/20
# content:    Test the package on same artificial data
import numpy as np
import pandas as pd
import anndata
import anndata_utils


def test_version():
    print(anndata_utils.version)
