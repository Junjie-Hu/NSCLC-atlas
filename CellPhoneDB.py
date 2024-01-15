import os
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

UMI = sc.read("all_sample_new.h5ad")
# get Endo
Endo = sc.read("Endo/Endo_cluster.h5ad")
#pick SELE+ Venule
SELE = Endo[Endo.obs.Subtype == "SELE+ Venule"]
# get macrophage
ann = pd.read_csv("MPs/MonoMacro_subtype_res.csv",index_col = 0)
# remove cycling cell
ann = ann.loc[ann.Subtype != "MKI67+ MoDM",:]
# neutrophil
neu = pd.read_csv("MPs/Neutrophil_subtype_res.csv",index_col=0)

# down sampling to 1/10
import random

ann.shape[0]/10
mpcb = random.sample(ann.index.tolist(), 24462)
ann_use = ann.loc[mpcb,:]
neu.shape[0]/10
neucb = random.sample(neu.index.tolist(), 1991)
neu_use = neu.loc[neucb,:]

Cell = SELE.obs.index.tolist() + ann_use.index.tolist() + neu_use.index.tolist()
Subtype = SELE.obs.Subtype.tolist() + ann_use.Subtype.tolist() + neu_use.Myeloid_lineage.tolist()
meta = pd.DataFrame(Cell,columns=['Cell'])
meta['cell_type'] = Subtype

meta.to_csv("cellphonedb/SELE_Myeloid_interction.csv",index = False)

sce = UMI[Cell,:]
sce.obs['cell_type'] = Subtype
sce.write("cellphonedb/SELE_Myeloid.h5ad")

# run cellphonedb v4
deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
        cpdb_file_path = "cellphonedb.zip",
        meta_file_path = "SELE_Myeloid_interction.csv",
        counts_file_path = "SELE_Myeloid.h5ad",
        counts_data = 'hgnc_symbol',
        output_path = "SELE_Myeloid",
        threshold = 0.1,
        threads = 10)

