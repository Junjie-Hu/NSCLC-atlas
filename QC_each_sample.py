import os

import pandas as pd 
import numpy as np
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt
from collections import Counter

### sample list
sams = os.listdir("./raw_matrix/matrix")
len(sams)

### QC each sample
i=0
sams[i]

## read matrix
mat = mmread("./raw_matrix/matrix/"+sams[i]+"/matrix.mtx")
gene = pd.read_csv("./raw_matrix/matrix/"+sams[i]+"/genes.tsv",sep="\t",header=None)
CB = pd.read_csv("./Singleron/raw_matrix/matrix/"+sams[i]+"/barcodes.tsv",sep="\t",header=None)
mat2 = mat.todense()
df = pd.DataFrame(mat2)
df["gene_symbols"] = gene[1].values
# sum UMIs for the genes with identical gene name
cols = CB.shape[0]
df_sum = df.groupby("gene_symbols", as_index=True)[list(df.columns[0:cols])].sum()
# add sample id to cell barcode
df_sum.columns = [ sams[i] +"_"+ a for a in CB[0].values]
Counter(gene[1].duplicated())
gene2 = gene.loc[[not i for i in gene[1].duplicated()],:]
gene2.index = gene2[1]
del gene2[1]
my_dict = gene2.to_dict()[0]
geneinfo = pd.DataFrame([my_dict[i] for i in df_sum.index],index=df_sum.index,columns=['genes_ids'])
# create scanpy object
sce = sc.AnnData(df_sum.T, var = geneinfo)

## scanpy processing
sce.var_names_make_unique()
sce.obs_names_make_unique()
# first set a cutoff value of 100 to get more neutrophils
sc.pp.filter_cells(sce, min_genes=100)
# mt genes
sce.var['mt'] = sce.var_names.str.match(r'MT-')
sc.pp.calculate_qc_metrics(sce, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# rb genes
sce.var['rb'] = sce.var_names.str.match(r'^RP[SL][0-9]')
sc.pp.calculate_qc_metrics(sce, qc_vars=['rb'], percent_top=None, log1p=False, inplace=True)
# remove cell with mt > 30 % or rb > 30%
sce = sce[sce.obs.pct_counts_mt <= 30, :]
sce = sce[sce.obs.pct_counts_rb <= 30, :]
# remove genes > 6000 and total counts > 30000
sce = sce[sce.obs.n_genes <= 6000,:]
sce = sce[sce.obs.total_counts <= 30000,:]

del sce.obs['total_counts_mt'],sce.obs['pct_counts_mt'],sce.obs['total_counts_rb'],sce.obs['pct_counts_rb']
del sce.var["mt"],sce.var["rb"]

## scrublet to remove doulets
scrub = scr.Scrublet(sce.to_df(), expected_doublet_rate=0.025)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
pd.value_counts(predicted_doublets)
sce.obs["doublet_scores"] = doublet_scores
sce.obs["predicted_doublets"] = predicted_doublets
# remove doublet predicted by scrublet
sce = sce[sce.obs["predicted_doublets"] == False,:]

## save
sce.write("../scanpy_raw/"+sams[i]+".h5ad")

## primary clustering to identify neutrophils
sce = sc.read("../scanpy_raw/"+sams[i]+".h5ad")
sc.pp.normalize_total(sce, target_sum=1e4) 
sc.pp.log1p(sce)
sc.pp.highly_variable_genes(sce, min_mean=0.0125, max_mean=3, min_disp=0.5)
sce.raw = sce
sce = sce[:, sce.var.highly_variable]
sc.pp.regress_out(sce, ['total_counts'])
sc.pp.scale(sce, max_value=10)
sc.tl.pca(sce, svd_solver='arpack')
sc.pl.pca_variance_ratio(sce, log=True)
sc.pp.neighbors(sce, n_neighbors=10, n_pcs=30)
sc.pp.scale(adata, max_value=10)
sc.tl.leiden(sce, resolution=0.8)
sc.pl.rank_genes_groups(sce, n_genes=15)
marker_genes = ["PTPRC","EPCAM",'SFTPC','SCGB3A1',"TPPP3","KRT17","CD3D","IL7R","CD8A","NKG7","CD79A","MS4A1","MZB1","JCHAIN",
                "LYZ","MARCO","CD14","CD68","FCN1",'VCAN','S100A8','FCGR3B',"KIT","MS4A2","DCN","TAGLN","RGS5","PECAM1","VWF",
               "PLD4","IL3RA","TOP2A","MKI67"]
sc.pl.dotplot(sce, marker_genes, groupby='leiden')
## annotate neutrophil
marker_genes = ["PTPRC","EPCAM",'SFTPC','SCGB3A1',"TPPP3","KRT17","CD3D","IL7R","CD8A","NKG7","CD79A","MS4A1","MZB1","JCHAIN",
                "LYZ","MARCO","CD14","CD68","FCN1",'VCAN','S100A8','FCGR3B',"KIT","MS4A2","DCN","TAGLN","RGS5","PECAM1","VWF",
               "PLD4","IL3RA","TOP2A","MKI67"]
sc.pl.dotplot(sce, marker_genes, groupby='leiden')

sce.obs['scanpy_lineage'] = "Other"
sce.obs.loc[sce[sce.obs.leiden.isin(["5"]),:].obs.index,"scanpy_lineage"] = "Neutrophil"
sce.obs.loc[sce[sce.obs.leiden.isin(["20"]),:].obs.index,"scanpy_lineage"] = "Doublets"
Counter(sce.obs.scanpy_lineage.values)
# qc
sce = sce[sce.obs.scanpy_lineage.values !="Doublets",:]
# neutrophil >=100 genes and other >= 300 genes
keep = []
for n in range(sce.obs.shape[0]):
    keep1 = sce.obs.n_genes[n] >= 300 or sce.obs.scanpy_lineage[n] == "Neutrophil"
    keep.append(keep1)
Counter(keep)
sce = sce[keep,:]

adata = sc.read("../scanpy_raw/"+sams[i]+".h5ad")
adata = adata[sce.obs.index,:]
## save
sce.write("../scanpy_raw/"+sams[i]+".h5ad")