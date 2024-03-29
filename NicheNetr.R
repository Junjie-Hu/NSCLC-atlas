library(Seurat)
library(ggplot2)
library(future)
library(tidyverse)
library(nichenetr)

RNA = readRDS("Endo.rds")

RNA.use <- subset(RNA,Subtype == "SELE+ Venule")
RNA_mat <- as.matrix(RNA.use@assays$RNA@data)

### Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks
ref_dir <- "/media/inspur/AS2150G2/HJJ/scrna/NicheNet/"
ligand_target_matrix <- readRDS(paste0(ref_dir,"ligand_target_matrix.rds"))
lr_network <- readRDS(paste0(ref_dir,"lr_network.rds"))
weighted_networks <- readRDS(paste0(ref_dir,"weighted_networks.rds"))

### Define the gene set of interest and a background of genes
M_marker <- read.csv("Endo_scanpy_leiden_res.csv",stringsAsFactors=F)
MP_marker <- M_marker %>%  filter(cluster== "SELE+ Venule") %>% filter(avg_logFC > 0.5)%>% dplyr::top_n(50, avg_logFC)
geneset_oi <- MP_marker$gene %>% .[. %in% rownames(ligand_target_matrix)]

### Define a background of genes
expressed_genes_receiver <- RNA_mat %>% apply(1,function(x){sum(x>0)/length(x)})  %>% .[. >= 0.1] %>% names()
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## define potential_ligands
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- ligands ## here we used all ligands intersect(ligands,expressed_genes_sender)
receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,expressed_genes_receiver)
lr_network_expressed <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
potential_ligands <- lr_network_expressed %>% pull(from) %>% unique()

### Perform NicheNet’s ligand activity analysis on the gene set of interest
ligand_activities <- predict_ligand_activities(geneset = geneset_oi, 
	background_expressed_genes = background_expressed_genes, 
	ligand_target_matrix = ligand_target_matrix, 
	potential_ligands = potential_ligands)
best_upstream_ligands <- ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

### Infer target genes of top-ranked ligands and visualize in a heatmap
active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique()
vis_ligand_target <- active_ligand_target_links[order_targets,order_ligands] %>% t()

MP_specific <-vis_ligand_target %>% make_heatmap_ggplot("Pro SELE+ Venule maturation ligands","SELE+ Venule-specific maturation signature", color = "red",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + 
scale_fill_gradient2(low = "whitesmoke",  high = "red", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

pdf("Endo_SELE_NicheNetr.pdf",width=7,height=6)
MP_specific
dev.off()

