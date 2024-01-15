library(Seurat)

RNA = readRDS("../Treg.rds")

# pick Treg
Treg = subset(RNA, Subtype != "CD8+ Treg")

st = read.csv("./stress_gene411.csv",header=T,stringsAsFactors = F)$gene

# Cytosig each subtype for each egimen
RNAlist = SplitObject(Treg, split.by = "batch")
sams = names(RNAlist)
for(i in sams){
    print(paste0("Runing ",i))
    RNA.tmp = RNAlist[[i]]
    data.mat = as.matrix(RNA.tmp@assays$RNA@data)
    data.sig = data.mat
    # 
    for(j in 1:nrow(data.mat)){
    data.sig[j,]=scale(data.mat[j,],center=T,scale=F)
    }
    ### remove 411 stress genens
    data.sig = data.sig[!c(rownames(data.sig) %in% st),]
    # SELL+ Treg
    SELL.num = which(RNA.tmp$Subtype == "SELL+ Treg")
    data.sig1 = data.sig[,SELL.num]
    SELL.run = length(SELL.num) >= 10
    if(SELL.run){write.csv(data.sig1,"Treg_SELL_a_sample_CytoSig.csv",quote=F)}
    # TNFRSF9+ Treg
    TNFRSF9.num = which(RNA.tmp$Subtype == "TNFRSF9+ Treg")
    data.sig2 = data.sig[,TNFRSF9.num]
    TNFRSF9.run = length(TNFRSF9.num) >= 10
    if(TNFRSF9.run){write.csv(data.sig2,"Treg_TNFRSF9_a_sample_CytoSig.csv",quote=F)}
    # NR4A2+ Treg
    NR4A2.num = which(RNA.tmp$Subtype == "NR4A2+ Treg")
    data.sig3 = data.sig[,NR4A2.num]
    NR4A2.run = length(NR4A2.num) >= 10
    if(NR4A2.run){write.csv(data.sig3,"Treg_NR4A2_a_sample_CytoSig.csv",quote=F)}
    ## cytosig
    cmd1 = "CytoSig_run.py -i Treg_SELL_a_sample_CytoSig.csv -o Treg_SELL_a_cytosig -s 1"
    if(SELL.run){system(cmd1)}
    cmd2 = "CytoSig_run.py -i Treg_TNFRSF9_a_sample_CytoSig.csv -o Treg_TNFRSF9_a_cytosig -s 1"
    if(TNFRSF9.run){system(cmd2)}
    cmd3 = "CytoSig_run.py -i Treg_NR4A2_a_sample_CytoSig.csv -o Treg_NR4A2_a_cytosig -s 1"
    if(NR4A2.run){system(cmd3)}
    ## calculate proliferation score for each cell
    cmd4 = paste0("python Cal_proliferation_score.py -i Treg_SELL_a_sample_CytoSig.csv -n False -o proliferate2/",i,"_SELL_prolifertion_score.csv")
    if(SELL.run){system(cmd4)}
    cmd5 = paste0("python Cal_proliferation_score.py -i Treg_TNFRSF9_a_sample_CytoSig.csv -n False -o proliferate2/",i,"_TNFRSF9_prolifertion_score.csv")
    if(TNFRSF9.run){system(cmd5)}
    cmd6 = paste0("python Cal_proliferation_score.py -i Treg_NR4A2_a_sample_CytoSig.csv -n False -o proliferate2/",i,"_NR4A2_prolifertion_score.csv")
    if(NR4A2.run){system(cmd6)}
    ## correlation for SELL Treg
    if(SELL.run){
    prof = read.csv(paste0("proliferate2/",i,"_SELL_prolifertion_score.csv"),header=T,stringsAsFactors = F)
    cytosig = read.table("Treg_SELL_a_cytosig.Zscore",header=T,stringsAsFactors = F,sep="\t",row.names=1)
    cytosig = as.data.frame(t(cytosig))
    res = list()
    for(k in colnames(cytosig)){
    a = cor.test(prof$Proliferation,cytosig[,k])
    b = data.frame(R = a$estimate,P = a$p.value)
    rownames(b) = k
    res[[k]] = b
    }
    res2 = do.call(rbind,res)
    res2$Q = p.adjust(res2$P,"fdr")
    write.csv(res2,paste0("prol_cytosig2/",i,"_SELL_res.csv"),quote=F)}
    ## correlation for TNFRSF9+ Treg
    if(TNFRSF9.run){
    prof = read.csv(paste0("proliferate2/",i,"_TNFRSF9_prolifertion_score.csv"),header=T,stringsAsFactors = F)
    cytosig = read.table("Treg_TNFRSF9_a_cytosig.Zscore",header=T,stringsAsFactors = F,sep="\t",row.names=1)
    cytosig = as.data.frame(t(cytosig))
    res = list()
    for(k in colnames(cytosig)){
    a = cor.test(prof$Proliferation,cytosig[,k])
    b = data.frame(R = a$estimate,P = a$p.value)
    rownames(b) = k
    res[[k]] = b
    }
    res2 = do.call(rbind,res)
    res2$Q = p.adjust(res2$P,"fdr")
    write.csv(res2,paste0("prol_cytosig2/",i,"_TNFRSF9_res.csv"),quote=F)}
    ## correlation for NR4A2+ Treg
    if(NR4A2.run){
    prof = read.csv(paste0("proliferate2/",i,"_NR4A2_prolifertion_score.csv"),header=T,stringsAsFactors = F)
    cytosig = read.table("Treg_NR4A2_a_cytosig.Zscore",header=T,stringsAsFactors = F,sep="\t",row.names=1)
    cytosig = as.data.frame(t(cytosig))
    res = list()
    for(k in colnames(cytosig)){
    a = cor.test(prof$Proliferation,cytosig[,k])
    b = data.frame(R = a$estimate,P = a$p.value)
    rownames(b) = k
    res[[k]] = b
    }
    res2 = do.call(rbind,res)
    res2$Q = p.adjust(res2$P,"fdr")
    write.csv(res2,paste0("prol_cytosig2/",i,"_NR4A2_res.csv"),quote=F)}
}
