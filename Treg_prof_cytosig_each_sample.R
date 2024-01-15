library(Seurat)

RNA = readRDS("../Treg.rds")

st = read.csv("./stress_gene411.csv",header=T,stringsAsFactors = F)$gene

# Cytosig each subtype for each egimen
RNAlist = SplitObject(Treg, split.by = "batch")
sams = names(RNAlist)
for(i in sams){
    print(paste0("Runing ",i))
    RNA.tmp = RNAlist[[i]]
    if(dim(RNA.tmp)[2] < 10){next}
    data.mat = as.matrix(RNA.tmp@assays$RNA@data)
    data.sig = data.mat
    # 
    for(j in 1:nrow(data.mat)){
    data.sig[j,]=scale(data.mat[j,],center=T,scale=F)
    }
    ### remove 411 stress genens
    data.sig = data.sig[!c(rownames(data.sig) %in% st),]
    write.csv(data.sig,"Treg_a_sample_CytoSig.csv",quote=F)
    ## cytosig
    cmd = "CytoSig_run.py -i Treg_a_sample_CytoSig.csv -o Treg_a_cytosig -s 1"
    system(cmd)
    ## calculate proliferation score for each cell
    cmd2 = paste0("python Cal_proliferation_score.py -i Treg_a_sample_CytoSig.csv -n False -o proliferate/",i,"_prolifertion_score.csv")
    system(cmd2)
    ## correlation
    prof = read.csv(paste0("proliferate/",i,"_prolifertion_score.csv"),header=T,stringsAsFactors = F)
    cytosig = read.table("Treg_a_cytosig.Zscore",header=T,stringsAsFactors = F,sep="\t",row.names=1)
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
    write.csv(res2,paste0("prol_cytosig/",i,"_res.csv"),quote=F)
}
