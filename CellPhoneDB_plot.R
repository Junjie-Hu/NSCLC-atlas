library(ktplots)
library(ggplot2)

setwd("./cellphonedb/SELE_Myeloid/")

pvalues <- read.delim("statistical_analysis_pvalues_09_13_2023_17:18:52.txt", check.names = FALSE)
means <- read.delim("statistical_analysis_means_09_13_2023_17:18:52.txt", check.names = FALSE)
sig.means <- read.delim("statistical_analysis_significant_means_09_13_2023_17:18:52.txt", check.names = FALSE)
deconvoluted <- read.delim("statistical_analysis_deconvoluted_09_13_2023_17:18:52.txt", check.names = FALSE)

plot_ord <- c("CD14+SELL+ Mono","CD14+CLEC10A+ Mono","CD14+CCR2+ Mono","CD16+CX3CR1+ Mono","FOLR2+ MoDM","CHIT1+ MoDM","CXCL3+ MoDM",
                "SPP1+ MoDM","CXCL9+ MoDM","Neutrophil")
cci_use1 = c("SELE+ Venule|CD14+SELL+ Mono","SELE+ Venule|CD14+CLEC10A+ Mono","SELE+ Venule|CD14+CCR2+ Mono","SELE+ Venule|CD16+CX3CR1+ Mono","SELE+ Venule|FOLR2+ MoDM","SELE+ Venule|CHIT1+ MoDM","SELE+ Venule|CXCL3+ MoDM",
               "SELE+ Venule|SPP1+ MoDM","SELE+ Venule|CXCL9+ MoDM","SELE+ Venule|Neutrophil")
cci_use2 <- c("CD14+SELL+ Mono|SELE+ Venule","CD14+CLEC10A+ Mono|SELE+ Venule","CD14+CCR2+ Mono|SELE+ Venule","CD16+CX3CR1+ Mono|SELE+ Venule","FOLR2+ MoDM|SELE+ Venule","CHIT1+ MoDM|SELE+ Venule","CXCL3+ MoDM|SELE+ Venule",
               "SPP1+ MoDM|SELE+ Venule","CXCL9+ MoDM|SELE+ Venule","Neutrophil|SELE+ Venule")
cci = c(cci_use1,cci_use2)

pvalues.use = pvalues[,c(colnames(pvalues)[1:11],cci)]
means.use = means[,c(colnames(sig.means)[1:11],cci)]

chemokines <- grep("^CXC|CX3C|CCL|CCR|CX3|XCL|XCR|IL[0-9]|CSF", means.use$interacting_pair,value = T)
adhesion <- grep("^VCAM|ICAM|SEL", means.use$interacting_pair,value = T)

pvalues.use2 = subset(pvalues.use, interacting_pair %in% c(chemokines,adhesion))
means.use2 = subset(means.use, interacting_pair %in% c(chemokines,adhesion))


dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = NULL,
                    scale_mean = FALSE,
                    width = 8,
                    height = 10,
                    means = all_means,
                    pvalues = all_pval,
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){

  #all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  #all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  all_pval = pvalues
  all_means = means
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }

  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  sel_means[sel_means == 0] = 1
  sel_means = log2(sel_means)
  ### scale means for each pair
  if(scale_mean){
    sel_means = as.data.frame(t(scale(t(sel_means))))
  }

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  #pr[pr==0] = 1
  #plot.data = cbind(plot.data,log2(pr))
  plot.data = cbind(plot.data,pr)
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  #my_palette <- colorRampPalette(c("blue", "green", "yellow", "red"), alpha=TRUE)(n=399)
  my_palette <- colorRampPalette(c("#228B22","YellowGreen","#FAFAD2", "#FFD700", "red"), alpha=TRUE)(n=399)

  ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn(ifelse(scale_mean,'Scaled Level','Level'), colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  if(!is.null(filename)){
     if (output_extension == '.pdf') {
      ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
      ggsave(filename, width = width, height = height, limitsize=F)
  }
  }
}

selected_rows = NULL
selected_columns = NULL
filename = "SELE_Myeloid"
scale_mean = FALSE
width = 8
height = 10
means = means.use2
pvalues = pvalues.use2
means_separator = '\t'
pvalues_separator = '\t'
output_extension = '.pdf'

 all_pval = pvalues.use2
  all_means = means.use2
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }

  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  sel_means[sel_means == 0] = 1
  sel_means = log2(sel_means)
  ### scale means for each pair
  if(scale_mean){
    sel_means = as.data.frame(t(scale(t(sel_means))))
  }

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  #pr[pr==0] = 1
  #plot.data = cbind(plot.data,log2(pr))
  plot.data = cbind(plot.data,pr)
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  my_palette <- colorRampPalette(c("#228B22","YellowGreen","#FAFAD2", "#FFD700", "red"), alpha=TRUE)(n=399)

  ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn(ifelse(scale_mean,'Scaled Level','Level'), colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

ggsave("SELE_Myeloid_plot.pdf", width = 8, height = 30, limitsize=F)


## select pairs
pair_plot = read.table("pair_plot.txt",header = F,stringsAsFactors = F)

pvalues.use3 = subset(pvalues.use2, interacting_pair %in% pair_plot$V1)
means.use3 = subset(means.use2, interacting_pair %in% pair_plot$V1)

selected_rows = NULL
selected_columns = NULL
filename = "SELE_Myeloid_select"
scale_mean = T
width = 8
height = 10
means = means.use3
pvalues = pvalues.use3
means_separator = '\t'
pvalues_separator = '\t'
output_extension = '.pdf'

all_pval = pvalues.use3
  all_means = means.use3
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }

sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

  if(scale_mean){
    sel_means = as.data.frame(t(scale(t(sel_means))))
    #sel_means = as.data.frame(scale(sel_means))
  }

sel_means[sel_means > 2] = 2

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  # pr[pr==0] = 1
  # plot.data = cbind(plot.data,log2(pr))
  plot.data = cbind(plot.data,pr)
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

plot.data$pair = factor(plot.data$pair,levels = rev(pair_plot$V1))

  my_palette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")), alpha=TRUE)(n=100)

  ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn(ifelse(scale_mean,'Scaled Level','Level'), colors=my_palette) +
  #scale_color_gradientn(ifelse(scale_mean,'Interaction \n Level','Level'), colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12, colour = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        ,axis.ticks.length = unit(0.2,"cm"),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

ggsave("SELE_Myeloid_select_plot.pdf", width = 11.5, height = 5.5, limitsize=F)
