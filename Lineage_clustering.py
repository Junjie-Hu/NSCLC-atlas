import os

import pandas as pd 
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from collections import Counter

sc._settings.ScanpyConfig.n_jobs = 10
sc.logging.print_header()
sc.set_figure_params()

## read h5ad
sce = sc.read("all_sample_new.h5ad")
sc.pl.highest_expr_genes(sce,n_top=20)
sce.var['mt'] = sce.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(sce, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(sce, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0, multi_panel=True)

## standard pipline
sc.pp.normalize_total(sce, target_sum=1e4) 
sc.pp.log1p(sce)

sc.pp.highly_variable_genes(sce, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key="batch")
Counter(sce.var.highly_variable)
sce.raw = sce
sce = sce[:, sce.var.highly_variable]
sc.pp.regress_out(sce, keys=['total_counts', 'pct_counts_mt'],n_jobs=30)
sc.pp.scale(sce, max_value=10)
sc.tl.pca(sce, svd_solver='arpack')
sc.pl.pca_variance_ratio(sce, log=True)
## BBKNN
sc.external.pp.bbknn(sce, batch_key='batch')
sc.pp.neighbors(sce, n_neighbors=10, n_pcs=40)
sc.tl.umap(sce)
sc.tl.leiden(sce, resolution=1)

marker_genes = ["PTPRC","EPCAM",'SFTPC',"TPPP3",'SCGB3A1',"KRT17","CD3D","IL7R","CD8A","NKG7","CD79A","MS4A1","MZB1","JCHAIN",
                "LYZ","MARCO","CD14","CD68","FCN1",'VCAN','S100A8','FCGR3B',"KIT","MS4A2","DCN","TAGLN","RGS5","PECAM1","VWF",
               "PLD4","IL3RA","TOP2A","MKI67"]
sc.pl.dotplot(sce, marker_genes, groupby='leiden')

## annotation lineage
# define fucntion
def listre(cluster,num,cell):
    for i in num:
        cluster[i] = cell
    return cluster

epi = [5,10,12,15,19,21,22,24,25,28,31,32,36,37,38,39,40,34]
Tc = [1,2,8,17,27]
NK = [18]
ILC = []
Bc = [0,33]
Pl = [3,20]
pDC = [35]
My = [4,6,9,13,16,30,42]
Ma = [11]
Ne = [26]
Fib = [14,23]
En = [7,41]
CI = [29,33]
Du = []
clus_num = len(set(sce.obs.leiden.values.tolist()))
free_ann = [None]*clus_num
free_ann = listre(free_ann,epi,"Epithelial")
free_ann = listre(free_ann,Tc,"T")
free_ann = listre(free_ann,NK,"NK")
free_ann = listre(free_ann,ILC,"ILC")
free_ann = listre(free_ann,Bc,"B")
free_ann = listre(free_ann,Pl,"Plasma")
free_ann = listre(free_ann,My,"Myeloid")
free_ann = listre(free_ann,Ma,"Mast")
free_ann = listre(free_ann,Ne,"Neutrophil")
free_ann = listre(free_ann,Fib,"Mesenchymal")
free_ann = listre(free_ann,En,"Endothelial")
free_ann = listre(free_ann,pDC,"pDC")
free_ann = listre(free_ann,CI,"Cycling immune")
free_ann = listre(free_ann,Du,"Doublets")

try:
	free_ann.index(None)
except Exception:
	print("all clusters have been annotated")

sce.obs['all_leiden'] = sce.obs.leiden.values
sce.obs['lineage'] = sce.obs.leiden.values 
new_cluster_names_unique=list(set(free_ann))
new_cluster_names_unique.sort(key=free_ann.index)
celltype = [free_ann[int(i)] for i in sce.obs.leiden.values.tolist()]
sce.obs.lineage = celltype
sce.obs.lineage = sce.obs.lineage.astype('category')
sce.obs.lineage.cat.set_categories(new_cluster_names_unique, inplace=True)

sce.write("all_cluster_new.h5ad")

## split cycling immune cell
cyc = sce[sce.obs.lineage == "Cycling immune",:]
sc.pp.neighbors(cyc, n_neighbors=10, n_pcs=20)
sc.tl.umap(cyc)
sc.tl.leiden(cyc, resolution=0.4)
# find marker
cyc.uns['log1p']['base'] = None
sc.tl.rank_genes_groups(cyc, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(cyc, n_genes=15)
# define
sce.obs.loc[cyc[cyc.obs.leiden.isin(["4","2","3"]),:].obs.index,"lineage"] = "Myeloid"
sce.obs.loc[cyc[cyc.obs.leiden.isin(["1","0"]),:].obs.index,"lineage"] = "T"
sce.obs.loc[cyc[cyc.obs.leiden.isin(["2"]),:].obs.index,"lineage"] = "B"
sce.obs.loc[cyc[cyc.obs.leiden.isin(["5"]),:].obs.index,"lineage"] = "Plasma"

sce.write("all_cluster_new.h5ad")

### split CD8+, CD4+, Treg lienage from all T cells
sce = sc.read("all_cluster_new.h5ad")
sce = sce.raw.to_adata()
# pick T/NK cells
sce = sce[sce.obe.lineage.isin(["T","NK"]),:]
# clustering
sce.uns['log1p']['base'] = None
sc.pp.highly_variable_genes(sce, min_mean=0.0125, max_mean=3, min_disp=0.5)
Counter(sce.var.highly_variable)
sce.raw = sce
sce = sce[:, sce.var.highly_variable]
sc.pp.regress_out(sce, keys=['total_counts'],n_jobs=30)
sc.pp.scale(sce, max_value=10)
sc.tl.pca(sce, svd_solver='arpack')
sc.pl.pca_variance_ratio(sce, log=True)
sc.external.pp.bbknn(sce, batch_key='batch')
sc.pp.neighbors(sce, n_neighbors=10, n_pcs=30)
sc.tl.umap(sce)
sc.tl.leiden(sce, resolution=1)
marker_genes = ["CD3D",'CD8A','CD8B','CD4','FOXP3','IL2RA',"NKG7","KLRD1","IL7R",'KIT']
sc.pl.dotplot(sce, marker_genes, groupby='leiden')
# annotation
sce.obs['TNK_lineage'] = sce.obs['lineage']
sce.obs.TNK_lineage.cat.set_categories(["CD8T","CD4T","Treg","NK","ILC3"], inplace=True)
sce.obs.loc[sce[sce.obs.leiden == "11",:].obs.index,"TNK_lineage"] = "ILC3"
sce.obs.loc[sce[sce.obs.leiden.isin(["5","6"]),:].obs.index,"TNK_lineage"] = "NK"
sce.obs.loc[sce[sce.obs.leiden.isin(["1"]),:].obs.index,"TNK_lineage"] = "Treg"
sce.obs.loc[sce[sce.obs.leiden.isin(["3","4"]),:].obs.index,"TNK_lineage"] = "CD8T"
sce.write("TNK_cluster_new.h5ad")

### split Monocyte/Macrophage and DC lienage from myeloid cells
sce = sc.read("all_cluster_new.h5ad")
sce = sce.raw.to_adata()
# pick myeloid cells
sce = sce[sce.obe.lineage.isin(["Myeloid","Neutrophil"]),:]
# clustering
sce.uns['log1p']['base'] = None
sc.pp.highly_variable_genes(sce, min_mean=0.0125, max_mean=3, min_disp=0.5)
Counter(sce.var.highly_variable)
sce.raw = sce
sce = sce[:, sce.var.highly_variable]
sc.pp.regress_out(sce, keys=['total_counts'],n_jobs=30)
sc.pp.scale(sce, max_value=10)
sc.tl.pca(sce, svd_solver='arpack')
sc.pl.pca_variance_ratio(sce, log=True)
sc.external.pp.harmony_integrate(sce,"batch")
sc.pp.neighbors(sce, n_neighbors=10, n_pcs=30)
sc.tl.umap(sce)
sc.tl.leiden(sce, resolution=0.8)
marker_genes = ["PTPRC","CD68","MERTK","CD14","FCGR3A","MRC1","CD1C",'LTB',"CCR7","CLEC9A","XCR1","GPNMB","TREM2","LAMP3",'SERPINB9',"CD2","BANK1",
               "FCN1","CD1A","CXCR4","CD83","CD86"]
sc.pl.dotplot(sce, marker_genes, groupby='leiden')
# annotation
sce.obs['Myeloid_lineage'] = sce.obs['lineage']
sce.obs.TNK_lineage.cat.set_categories(["Mono/Macro","DC","Neutrophil","Doublets"], inplace=True)
sce.obs.loc[sce[sce.obs.leiden.isin(["0","1","2","4","5","7","8","10"]),:].obs.index,"Myeloid_lineage"] = "Mono/Macro"
sce.obs.loc[sce[sce.obs.leiden.isin(["3","11","14"]),:].obs.index,"Myeloid_lineage"] = "DC"
sce.obs.loc[sce[sce.obs.leiden.isin(["7"]),:].obs.index,"Myeloid_lineage"] = "Neutrophil"
sce.obs.loc[sce[sce.obs.leiden.isin(["6","9","12","13"]),:].obs.index,"Myeloid_lineage"] = "Doublets"
sce.write("Myeloid_cluster_new.h5ad")

### re-annotation T and myeloid cells
sce = sc.read("all_cluster_new.h5ad")
TNK = sc.read("TNK_cluster_new.h5ad")
Myeloid = sc.read("Myeloid_cluster_new.h5ad")
new_celltype = list(set(sce.obs.lineage)) + ["Mono/Macro","DC","CD8T","CD4T","Treg","ILC3"]
sce.obs.lineage.cat.set_categories(new_celltype, inplace=True)

sce.obs.loc[TNK.obs.index,"lineage"] = TNK.obs["TNK_lineage"]
sce.obs.loc[Myeloid.obs.index,"lineage"] = Myeloid.obs["Myeloid_lineage"]

sce = sce[sce.obs.lineage != "Doublets",:]
new_celltype2 = list(set(sce.obs.lineage))
sce.obs.lineage.cat.set_categories(new_celltype2, inplace=True)

sce.write("all_cluster_new.h5ad")

