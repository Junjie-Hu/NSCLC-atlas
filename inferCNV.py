import os

import pandas as pd 
import numpy as np
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
from collections import Counter

sc.settings.set_figure_params(figsize=(5, 5))
sc.logging.print_header()

# pick non-immune cells
sceall = sc.read("all_cluster_new.h5ad")
sceall = aceall.raw.to_adata()
sceall = sceall[sceall.obs.lineage.isin("Epithelial","Mesenchymal","Endothelial"),:]

cnv.io.genomic_position_from_gtf(gtf_file="../index/gencode.v38.annotation.gtf",
                                 adata=sceall,
                                 gtf_gene_id="gene_name")

# perform infercnv for each sample
i=0
sam = sceall.obs.batch.cat.categories[i]
sce = sceall[sceall.obs.batch == sam,:]

cnv.tl.infercnv(
   sce,
   reference_key="Subtype",
   reference_cat=["Mesenchymal","Endothelial"],
   window_size=250,
)
cnv.pl.chromosome_heatmap(sce, groupby="Subtype")
cnv.tl.pca(sce)
cnv.pp.neighbors(sce)
cnv.tl.leiden(sce)
cnv.tl.cnv_score(sce)
cnv.pl.chromosome_heatmap(sce, groupby="cnv_leiden", dendrogram=True,save=sam+"_cnvleiden_heat")
## boxplot
plt.figure(figsize=(10,5))
b=[sce.obs.cnv_score[sce.obs.cnv_leiden == str(i)] for i in range(0,13)]
plt.boxplot(b)

marker_genes = ["EPCAM","NAPSA","NKX2-1","KRT7","MUC1",'SFTPC',"SFTPA1","PGC","TPPP3",'SCGB3A1',"KRT6A","KRT5","TP63","SOX2",
                "DCN","TAGLN","RGS5","PECAM1","VWF","GRP","CHGA","SYP","NCAM1","ASCL1","NEUROD1","SST","VIM","GPX2","MKI67"]
sc.pl.dotplot(sce, marker_genes, groupby='cnv_leiden')
#Classifying tumor cells
sce.obs["cnv_status"] = "Normal"
sce.obs.loc[sce.obs["cnv_leiden"].isin(["5","4"]),"cnv_status"] = "Malignant"
sce[sce.obs.Subtype == "Epithelial",:].obs[["batch","sample_lineage","lineage","cnv_leiden","leiden","cnv_score","cnv_status"]].to_csv("./infercnv/"+sam+"_infercnvpy_res.csv")


