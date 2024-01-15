import os

import pandas as pd 
import numpy as np
import scanpy as sc

## read metadata
meta = pd.read_csv("metadata.csv")

sams = meta["Sample"].tolist()
cans = meta["Sample"].tolist()
len(sams);len(cans)

for i in range(0,len(sams)):
    ### raw scanpy
    adata = sc.read("scanpy_raw/"+sams[i]+".h5ad")
    del adata.obs["n_genes_by_counts"],adata.obs["total_counts"],adata.obs["doublet_scores"],adata.obs["predicted_doublets"]
    del adata.var['n_cells_by_counts'],adata.var['mean_counts'],adata.var['pct_dropout_by_counts'],adata.var['total_counts']
    cans[i] = adata

## merge
sce_all = cans[0].concatenate(cans[1:],join='outer',index_unique=None,batch_categories=sams)
# clean anndata.var
ENSG = []
for i in sce_all.var.index:
    a = np.array(list(set(sce_all.var.loc[i,:].values))).tolist()
    a.sort()
    ENSG.append(a[0])
sce_all.var["gene_ids"] = ENSG
sce_all.var = sce_all.var.drop(sce_all.var.columns.values.tolist()[0:120],axis=1)

sce_all.obs_names_make_unique()
sce_all.var_names_make_unique()


## remove genes that expressed in less than 30 cells
sc.pp.filter_genes(sce_all, min_cells=30)

## save
sce_all.write("all_sample_new.h5ad")