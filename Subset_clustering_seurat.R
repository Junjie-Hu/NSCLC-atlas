library(Seurat)
library(ggplot2)
library(dplyr)
library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 100*1024^3)

RNA = readRDS("CD8_new.rds")

## remove TCR genes regarding variable region and BCR genes
genes = rownames(RNA)
TCR = grep("^TR[ABDG][VDJ]",genes)
BCR = grep("^IG[HKLJ]",genes)

RNA = RNA[-c(TCR,BCR),]

RNA <- NormalizeData(object = RNA)

aRNA = RNA
# remove cycling leiden cluster to get HVG
aRNA = subset(aRNA, immune_leiden != 11)
aRNA = subset(aRNA, leiden != 13)
## use chr1-22 and X gene
noY = read.table("/media/inspur/AS2150G2/HJJ/scrna/chr1_22_X_genename.txt",header=T,stringsAsFactors = F)
aRNA = aRNA[rownames(aRNA) %in% noY$gene_name,]

aRNA <- FindVariableFeatures(object = aRNA,nfeatures = 2000)
# remove bad genes from HVGs
gene_ignore = read.csv("gene_ignore_T.csv",header=T,stringsAsFactors=F)
bkgene = c(gene_ignore$gene, "MALAT1","XIST","NEAT1")
hvg = aRNA@assays$RNA@var.features[!c(aRNA@assays$RNA@var.features %in% bkgene)]
RNA@assays$RNA@var.features <- hvg

RNA <- ScaleData(object = RNA,vars.to.regress=c("nCount_RNA"),features = RNA@assays$RNA@var.features)
RNA <-  RunPCA(RNA, features = RNA@assays$RNA@var.features)
library(harmony)
RNA <- RunHarmony(RNA,"batch",plot_convergence = F)
RNA <- FindNeighbors(object = RNA,  reduction = "harmony",verbose = T, dims = 1:20)
RNA <- FindClusters(object = RNA, verbose = T, resolution = 0.8)
RNA <- RunUMAP(RNA, reduction = "harmony", dims = 1:20)

# Find markers for each cluster
markers <- FindAllMarkers(object = RNA, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers, file=paste0("CD8_harmony_res0.8.csv"),quote=F)

# remove doublets 
RNAharm = subset(RNAharm,RNA_snn_res.0.8 != 12) # TPSB2

# re-do clustering above

## annoation
free_annotation = c("GZMK+ Tem","ZNF683+ Trm","NR4A1+ Tem","CD8+TCF7+ Tn","IFNG+ Tem","CXCL13+ Tex","XCL1+ Trm",
                   "KLRB1+ Tem","MKI67+ Tem","SYNE2+ Tem","TRDC+ gdT")
RNA@meta.data$Subtype <- plyr::mapvalues(x = RNA@meta.data$RNA_snn_res.0.8, from = 0:10, to = free_annotation)

saveRDS(RNA,"CD8_harmony_res0.8.rds")

### plot heatmap of marker genes
source("./Heat_Dot_data.R")
library(tidyverse)

genes_to_check = c("TRDC","TRGC1","SELL","LEF1","TCF7","CCR7","IL7R",
                   "ZNF683","ITGAE","GZMA",'GZMB',"GZMK",'IFNG',"PRF1",'GNLY',"NKG7","LAG3","TIGIT","PDCD1","HAVCR2","LAYN","ENTPD1","CTLA4",
                   "LTB","IL22","CXCL13","XCL1","XCL2","CCL3","CCL4","CCL5","CCL20","CD69","CD160","LAIR2","TNFRSF4","TNFRSF9","KLRB1","CXCR4","CXCR6",
                  "EOMES","MAF","TOX2","NR4A1","NR4A2","NR4A3","BATF",
                   "MKI67","TOP2A")
row_split = c(rep("gdT",2),rep("Naive",5),rep("TRM",2),rep("Cytotoxic",7),rep("Exhausted",7),rep("Chemokines",9),rep("Receptors",8),rep("TFs",7),rep("Other",2))
row_split = factor(row_split,levels = c("gdT","Naive","TRM","Cytotoxic","Exhausted","Chemokines","Receptors","TFs","Other"))
data.plot <- Heat_Dot_data(object=RNA,features=genes_to_check,group.by="Subtype")
plot_ord <- c("CD8+TCF7+ Tn","GZMK+ Tem","NR4A1+ Tem","IFNG+ Tem","KLRB1+ Tem","ZNF683+ Trm","XCL1+ Trm",
                   "SYNE2+ Tem","MKI67+ Tem","CXCL13+ Tex","TRDC+ gdT")
exp.mat <- data.plot %>% select(features.plot,id,avg.exp.scaled) %>% spread(id,avg.exp.scaled)
rownames(exp.mat) <- exp.mat$features.plot
exp.mat$features.plot <- NULL
exp.mat <- exp.mat[,plot_ord]
per.mat <- data.plot %>% select(features.plot,id,pct.exp) %>% spread(id,pct.exp)
rownames(per.mat) <- per.mat$features.plot
per.mat$features.plot <- NULL
per.mat <- per.mat[,plot_ord]/100

max(exp.mat);min(exp.mat)

library(ComplexHeatmap)
library(circlize) ## color 
col_fun <- colorRamp2(c(-1.5, 0, 2.5), c("#0077b6", "#FFFFCC", "#e63946"))
# split heatmap
col_split = c(rep("ab",10),"gd")
col_split =factor(col_split,levels = c("ab","gd"))
# left annotation
c("gdT","Naive","TRM","Cytotoxic","Exhausted","Chemokines","Receptors","TFs","Other")
ha = HeatmapAnnotation(df = data.frame(Marker=row_split),which = "row",
                       col = list(Marker = c("gdT"= "#e76f51","Naive" = "#ffafcc","TRM"="#264653","Cytotoxic"="#0077b6","Exhausted"="#ddbea9",
                                          "Chemokines"="#57cc99","Receptors"="#dc2f02","TFs"="#fca311","Other"="#ff34a1")))
# top boxplot of n_genes
boxlist = list()
for(i in plot_ord){
    boxlist[[i]] = subset(RNA@meta.data,Subtype == i)$nFeature_RNA
}
ha_top = HeatmapAnnotation(nFeature_RNA = anno_boxplot(boxlist, which ="column",height = unit(2, "cm"),
                                         gp=gpar(fill = cols), outline = FALSE))

pdf("CD8T_marker_heatmap.pdf",width = 5.8,height = 16)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(col="grey"),
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,column_split = col_split,
        row_gap = unit(2, "mm"),column_gap =unit(2, "mm"), 
        left_annotation = ha,top_annotation = ha_top,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()

pdf("CD8T_marker_heatmap_dot.pdf",width = 5.8,height = 16)
Heatmap(exp.mat, col = col_fun,cluster_columns = F,cluster_rows = F,
        show_column_names = T,show_row_names = T,rect_gp=gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fil){
          grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y,r=per.mat[i,j]/2 * max(unit.c(width, height)),gp = gpar(fill = col_fun(exp.mat[i, j]), col = NA))},
        column_names_side = "top",row_names_side = "right",
        row_split = row_split,column_split = col_split,
        row_gap = unit(2, "mm"),column_gap =unit(2, "mm"), 
        left_annotation = ha,top_annotation = ha_top,
        heatmap_legend_param=list(title = "Expression",legend_height=unit(3, "cm")))
dev.off()
