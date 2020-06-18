[![Build Status](https://travis-ci.org/iosonofabio/anndata_utils.svg?branch=master)](https://travis-ci.org/iosonofabio/anndata_utils)

# anndata_utils
A bunch of odds and evens I wish anndata had.

```python
import anndata
import anndata_utils

adata = anndata.read_loom('dataset1.loom')

# Split into separate AnnData objects by metadata
adata_dic = anndata_utils.split(adata, 'obs', 'cellType')
```

Have fun!

