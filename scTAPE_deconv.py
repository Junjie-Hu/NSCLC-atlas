import re

import pandas as pd
import scanpy as sc
from TAPE import Deconvolution


sc_ref = pd.read_csv("All_celltype_scTAPE_count.txt",sep = "\t")

# TCGA bulk pre-process
ref = "./index/output.annotation"
anno = pd.read_csv(ref, delimiter='\t',header=None,names=["transcript_id","ensembl_id","gene_name"])
anno['ensembl_id'] = [re.sub("\\.[0-9]*$","",x)  for x in anno['ensembl_id']]

def Ensembl2name(df,colname,anno):
    annot_dict = anno.set_index('ensembl_id').to_dict()['gene_name']
    df[colname] = df[colname].apply(lambda x: annot_dict[x] if x in annot_dict.keys() else x)
    cols = df.shape[1]
    df_redup = df.groupby(colname, as_index=False)[list(df.columns[1:cols])].sum() # sum or mean
    return df_redup

bulk = pd.read_csv("TCGA-LUAD.htseq_fpkm.tsv",sep="\t")
bulk['Ensembl_ID'] = [re.sub("\\.[0-9]*$","",x)  for x in bulk['Ensembl_ID']]

bulk_redup = Ensembl2name(bulk, "Ensembl_ID", anno)
bulk_redup.index = bulk_redup["Ensembl_ID"]
del bulk_redup["Ensembl_ID"]
bulk_redup = bulk_redup.T
bulk_redup.to_csv("TCGA_LUAD_FPKM_trans.txt",sep="\t")

## run scTAPE
sc_ref = "./All_celltype_scTAPE_count.txt"
bulkdata = "./TCGA_LUAD_FPKM_trans.txt"

SignatureMatrix, CellFractionPrediction = \
    Deconvolution(sc_ref, bulkdata, sep='\t', scaler='mms',
                  datatype='counts', genelenfile='./GeneLength.txt',
                  mode='overall', adaptive=True, variance_threshold=0.98,
                  save_model_name=None,
                  batch_size=128, epochs=128, seed=1)

SignatureMatrix.to_csv("All_celltype_split_CAF_signature.csv")
CellFractionPrediction.to_csv("TCGA_LUAD_all_celltype_split_CAF_frac.csv")

