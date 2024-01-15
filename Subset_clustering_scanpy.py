import os

import pandas as pd 
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from collections import Counter

sc.logging.print_header()
sc.set_figure_params()
sc.set_figure_params(color_map = "viridis_r")

sce = sc.read("Endo_cluster_new.h5ad")
sce = sce.raw.to_adata()

asce = sce.copy()
## use chr1-22 and X gene
noY = pd.read_csv("/media/inspur/AS2150G2/HJJ/scrna/chr1_22_X_genename.txt",sep="\t")
gene_use = [i in noY["gene_name"].values for i in asce.var.index]
Counter(gene_use)
####
asce.uns['log1p']['base'] = None
sc.pp.highly_variable_genes(asce, min_mean=0.0125, max_mean=3, min_disp=0.5)

### remove bad genes from HVGs
bkgene = pd.read_csv("./gene_ignore_T.csv")
bklist = bkgene['gene'].values.tolist() + ["MALAT1","XIST","NEAT1"]
len(bklist)
asce.var.loc[asce.var.index.isin(bklist),"highly_variable"] = False
asce.raw = asce
asce = asce[:, asce.var.highly_variable]

sce.raw = sce
sce = sce[:, asce.var.index]

sc.pp.regress_out(sce, ['total_counts'],n_jobs = 30)
sc.pp.scale(sce, max_value=10)
sc.tl.pca(sce, svd_solver='arpack')
sc.pl.pca_variance_ratio(sce, log=True)
# harmony
sc.external.pp.harmony_integrate(sce,"batch")
sc.pp.neighbors(sce, n_neighbors=10, n_pcs=30,use_rep = "X_pca_harmony")
sc.tl.umap(sce)
sc.tl.leiden(sce, resolution=0.6)
# find marker
sce.uns['log1p']['base'] = None
sc.tl.rank_genes_groups(sce, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(sce, n_genes=15)

endo_genes = ["ITGB1","PECAM1","GJA5","DKK2","FBLN5","ACKR1","CPE","FENDRR","TBX3","CA4","CD36","FCN3","EDNRB","IL1RL1","APLN","TBX2",
                   "PROX1","TFF3","CCL21","PDPN","CXCR4","PGF","INSR","ESM1","CHST4", "FUT7", "LIPG", "ENPP2", "LIFR", "CH25H", "C4BPA", "ENPP6",
                   "SPARC","COL4A1","COL15A1","CD34","LAMA4","PXDN","LOX","SELE","SELP","VCAM1","ICAM1","CCL2","IL6","CSF3","IL32",
                   "CXCL9",'CXCL10','CXCL11',"KDR","C7","APLNR","CLU","OLFM1","CCL23"]
sc.pl.dotplot(sce,endo_genes, groupby='leiden')

sce.write("Endo_cluster_new.h5ad")

# remove doublets
sce = sce[sce.obs.leiden != "6",:]# LYZ, TROPBP, PTPRC
## clustering again

#### annotation
sce.obs['Subtype'] = "CPE- Venule"

sce.obs.loc[sce.obs.leiden.isin(["1"]),'Subtype'] = "CPE+ Venule"
sce.obs.loc[sce.obs.leiden.isin(["2"]),'Subtype'] = "PGF- Tip"
sce.obs.loc[sce.obs.leiden.isin(["3"]),'Subtype'] = "SELE+ Venule"
sce.obs.loc[sce.obs.leiden.isin(["9"]),'Subtype'] = "APLNR+ Venule"
sce.obs.loc[sce.obs.leiden.isin(["5"]),'Subtype'] = "Artery"
sce.obs.loc[sce.obs.leiden.isin(["4"]),'Subtype'] = "gCap"
sce.obs.loc[sce.obs.leiden.isin(["6"]),'Subtype'] = "PGF+ Tip"
sce.obs.loc[sce.obs.leiden.isin(["7"]),'Subtype'] = "Lymphatic EC"
sce.obs.loc[sce.obs.leiden.isin(["8"]),'Subtype'] = "Aerocyte"
sce.obs.loc[sce.obs.leiden.isin(["10"]),'Subtype'] = "MKI67+ EC"
sce.write("Endo_cluster_new.h5ad")

# split cycling  cell
cyc = sce[sce.obs.leiden == "10",:]
sc.pp.neighbors(cyc, n_neighbors=10, n_pcs=20)
sc.tl.umap(cyc)
sc.tl.leiden(cyc, resolution=0.5)
# annotation
sce.obs.loc[cyc[cyc.obs.leiden.isin(["1","5","0"]),:].obs.index,"Subtype"] = "CPE- Venule"
sce.obs.loc[cyc[cyc.obs.leiden.isin(["3"]),:].obs.index,"Subtype"] = "gCap"
sce.obs.loc[cyc[cyc.obs.leiden.isin(["4","2"]),:].obs.index,"Subtype"] = "PGF- Tip"

subtype = list(set(sce.obs.Subtype))
sce.obs.Subtype.cat.set_categories(subtype, inplace=True)

sce.write("Endo_cluster_new.h5ad")

## save metadata and UMAP embedding for marker heatmap in R
sce.obs.loc[:,["batch",'sample_lineage','lineage',"leiden","Subtype"]].to_csv("/Endo_scanpy_subtype_res.csv")
umap = pd.DataFrame(sce.obsm["X_umap"],columns = ['UMAP_1',"UMAP_2"])
umap.index = sce.obs.index
umap.to_csv("./Endo_scanpy_umap.csv")

