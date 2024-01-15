library(Seurat)
library(ggplot2)
library(dplyr)

RNA = readRDS("Mesen.rds")

# Cytosig for each subtype 
RNAlist = SplitObject(RNA, split.by = "Subtype")
sams = names(RNAlist)
resall = list()
for(i in sams){
    RNA.tmp = RNAlist[[i]]
    data.exp = as.matrix(RNA.tmp@assays$RNA@data)
    resall[[i]]=apply(data.exp,1,mean)
}
data.mat <- do.call(what = cbind,args=resall)

data.sig = data.mat
for(i in 1:nrow(data.mat)){
    data.sig[i,]=scale(data.mat[i,],center=T,scale=F)
}

# remove 411 stress genens
st = read.csv("./stress_gene411.csv",header=T,stringsAsFactors = F)$gene
data.sig = data.sig[!c(rownames(data.sig) %in% st),]
dim(data.sig)
write.csv(data.sig,"Mesen_subtype_CytoSig.csv",quote=F)

## bash 
# run cytosig
cmd = "CytoSig_run.py -i Mesen_subtype_CytoSig.csv -o Mesen_cytosig -s 1" 
system(cmd)

### CytoSig heatmap
score = read.table("Mesen_cytosig.Zscore",sep="\t",check.names = F)

library(ComplexHeatmap)
library(circlize) ## color
library(RColorBrewer)

cols = c(brewer.pal(12,"Set3"))
names(cols) = c("MYH11-PTN+ Pericyte","AA","MYH11-PTN- Pericyte","FAP+aSMA+ CAF","MYH11+ Pericyte","BCHE- SMC","ADH1B+ CAF","FAP+aSMA- CAF","CCL21+ PvC",
             "BCHE+ SMC","Mesothelial-like fib","CLU+ fib")
cols
annotation_colors = list(Subtype=cols)
ha = HeatmapAnnotation(df = data.frame(Subtype = prod),which = "row",
                       col = list(Subtype=cols))
col_fun2 <- colorRamp2(c(-16, -8,0, 8,16), c("#483D8B","#0000FF", "#FFFFF0", "red","darkred"))
ht = Heatmap(t(score2),col = col_fun,show_column_names = T,show_row_names = T,
        cluster_columns = T,cluster_rows = T,
        heatmap_legend_param=list(title = "Signaling Activity",legend_height=unit(3, "cm")),
        row_names_side = "left",column_dend_side = "bottom",column_names_side = "top",
             left_annotation = ha)
pdf("Mesen_chemokine_cytosig.pdf",width = 15,height = 4.2) 
draw(ht,  heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()

#### Cytosig for each subtype in each sample
RNA@meta.data$Group = paste(RNA$batch,RNA$Subtype,sep="_")
RNAlist = SplitObject(RNA, split.by = "Group")
sams = names(RNAlist)
resall = list()
for(i in sams){
    RNA.tmp = RNAlist[[i]]
    data.exp = as.matrix(RNA.tmp@assays$RNA@data)
    resall[[i]]=apply(data.exp,1,mean)
}
data.mat <- do.call(what = cbind,args=resall)
data.sig = data.mat
for(i in 1:nrow(data.mat)){
    data.sig[i,]=scale(data.mat[i,],center=T,scale=F)
}
# remove 411 stress genens
st = read.csv("./stress_gene411.csv",header=T,stringsAsFactors = F)$gene
data.sig = data.sig[!c(rownames(data.sig) %in% st),]
dim(data.sig)
write.csv(data.sig,"Mesen_subtype_each_sample_CytoSig.csv",quote=F)

## bash
# run cytosig
cmd = "CytoSig_run.py -i Mesen_subtype_each_sample_CytoSig.csv -o Mesen_each_sample_cytosig -s 1" 
system(cmd)

# CytoSig res
score = read.table("Mesen_each_sample_cytosig.Zscore",sep="\t",check.names = F)
score2 = t(score) %>% as.data.frame()
df = as.data.frame(table(RNA@meta.data$Group))
colnames(df) = c("Group","CellNumer")
rownames(df) = df$Group
df$Sample = gsub("_.*$","",df$Group)
df$Subtype = gsub("^P.*_","",df$Group)
df$Group = NULL
df = df[,c(2,3,1)]
head(df)
score3 = merge(score2,df, by=0)
score3 = score3[,c(1,53:55,2:52)]
head(score3);dim(score3)
write.csv(score3[,-1],"Mesen_cytokine_activity_each_subtype_each_sample.csv",quote=F,row.names=F)

