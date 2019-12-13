# scRNA Covariance Patient Classification
# 11/4/19

library(dplyr)
library(tibble)
library(Seurat)
library(corrplot)
library(caret)

# # load in data
# load("~/rao_lab/human_interactome/Final_Tissue_Merge_Batch_Fil.RData")
# 
# # All genes in object
# all_genes <- rownames(Tissue_Merge)
# 
# # Find common cell types for patients (Simple_Relabelling)
# # Don't use 181658
# # Use Myeloid, Epithelial, Fibroblast, T Cell, NK
# table(Tissue_Merge$Simple_Relabelling, Tissue_Merge$Patient_ID) > 0
# wanted_cell_type <- c('Myeloid', 'Epithelial', 'Fibroblast', 'T Cell', 'NK')
# 
# # Subset wanted cell types into new object
# Idents(Tissue_Merge) <- 'Simple_Relabelling'
# tissue_merge_class <- subset(Tissue_Merge, idents = wanted_cell_type)
# 
# # Run DE with batch correction between groups for each cell type 
# tissue_merge_class_DE <- FindAllMarkers(tissue_merge_class, test.use = 'LR', latent.vars = 'Patient_ID')
# 
# # Filter DE genes of RP and MT- genes
# filt_genes <- unique(tissue_merge_class_DE$gene[grep('^RP|^MT-', DE_filt$gene)])
# DE_filt <- tissue_merge_class_DE %>% group_by(cluster) %>% filter(!gene %in% filt_genes)
# 
# # Get genes where p_val_adj < .05
# DE_filt <- DE_filt %>% group_by(cluster) %>% filter(p_val_adj < .05)
# cell_type_gene_list <- lapply(wanted_cell_type, function(x) DE_filt$gene[which(DE_filt$cluster == x)])
# names(cell_type_gene_list) <- wanted_cell_type
# 
# # Create data frame of DE genes and wanted cell types
# tissue_data <- FetchData(Tissue_Merge, vars = c(unique(DE_filt$gene), 'Simple_Relabelling', 'Patient_ID', 'State'))
# tissue_data <- tissue_data %>% filter(Simple_Relabelling %in% wanted_cell_type & Patient_ID != '181658')
# 
# # Myeloid (> .8)
# myeloid_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[1]) %>% select(cell_type_gene_list[[1]], 'Simple_Relabelling', 'Patient_ID', 'State')
# myeloid_cor <- cor(myeloid_data[,1:length(cell_type_gene_list[[1]])], method = 'spearman')
# table(abs(myeloid_cor) > .7 & abs(myeloid_cor) != 1)
# corrplot(abs(myeloid_cor[1:30, 1:30]) > .7 & abs(myeloid_cor[1:30, 1:30]) < 1)
# 
# myeloid_cor_logic <- myeloid_cor > .7 & myeloid_cor < 1
# valid_n <- sum(myeloid_cor_logic)
# myeloid_cor_1 <- data.frame(row_gene = rep(0, times = valid_n), col_gene = rep(0, times = valid_n))
# n <- 1
# 
# for (i in 1:nrow(myeloid_cor)) {
#   for(j in 1:ncol(myeloid_cor)) {
#     if (myeloid_cor_logic[i,j]) {
#       myeloid_cor_1$row_gene[n] <- rownames(myeloid_cor_logic)[i]
#       myeloid_cor_1$col_gene[n] <- colnames(myeloid_cor_logic)[j]
#       n <- n + 1
#     }
#   }
# }
# 
# myeloid_genes <- unique(myeloid_cor_1$row_gene)
# 
# # Epithelial (> .7)
# epithelial_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[2]) %>% select(cell_type_gene_list[[2]], 'Simple_Relabelling', 'Patient_ID', 'State')
# epithelial_cor <- cor(epithelial_data[,1:length(cell_type_gene_list[[2]])], method = 'spearman')
# table(abs(epithelial_cor) > .7 & abs(epithelial_cor) != 1)
# corrplot(abs(epithelial_cor[1:30, 1:30]) > .7 & abs(epithelial_cor[1:30, 1:30]) < 1)
# 
# epithelial_cor_logic <- epithelial_cor > .7 & epithelial_cor < 1
# valid_n <- sum(epithelial_cor_logic)
# epithelial_cor_1 <- data.frame(row_gene = rep(0, times = valid_n), col_gene = rep(0, times = valid_n))
# n <- 1
# 
# for (i in 1:nrow(epithelial_cor)) {
#   for(j in 1:ncol(epithelial_cor)) {
#     if (epithelial_cor_logic[i,j]) {
#       epithelial_cor_1$row_gene[n] <- rownames(epithelial_cor_logic)[i]
#       epithelial_cor_1$col_gene[n] <- colnames(epithelial_cor_logic)[j]
#       n <- n + 1
#     }
#   }
# }
# 
# epithelial_genes <- unique(epithelial_cor_1$row_gene)
# 
# # Fibroblast (> .9)
# fibroblast_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[3]) %>% select(cell_type_gene_list[[3]], 'Simple_Relabelling', 'Patient_ID', 'State')
# 
# # Check gene standard deviation is > 0
# fibroblast_sd <- apply(fibroblast_data[,1:length(cell_type_gene_list[[3]])], 2, sd)
# x <- fibroblast_sd > 0
# fibroblast_data <- fibroblast_data[,x]
# 
# fibroblast_cor <- cor(fibroblast_data[,1:sum(x)], method = 'spearman')
# 
# fibroblast_cor_logic <- fibroblast_cor > .7 & fibroblast_cor < 1
# valid_n <- sum(fibroblast_cor_logic)
# fibroblast_cor_1 <- data.frame(row_gene = rep(0, times = valid_n), col_gene = rep(0, times = valid_n))
# n <- 1
# 
# for (i in 1:nrow(fibroblast_cor)) {
#   for(j in 1:ncol(fibroblast_cor)) {
#     if (fibroblast_cor_logic[i,j] & !is.na(fibroblast_cor_logic[i,j])) {
#       fibroblast_cor_1$row_gene[n] <- rownames(fibroblast_cor_logic)[i]
#       fibroblast_cor_1$col_gene[n] <- colnames(fibroblast_cor_logic)[j]
#       n <- n + 1
#     }
#   }
# }
# 
# fibroblast_genes <- unique(fibroblast_cor_1$row_gene)
# 
# # T Cell (> .8)
# t_cell_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[4]) %>% select(cell_type_gene_list[[4]], 'Simple_Relabelling', 'Patient_ID', 'State')
# 
# # Check gene standard deviation is > 0
# t_cell_sd <- apply(t_cell_data[,1:length(cell_type_gene_list[[4]])], 2, sd)
# x <- t_cell_sd > 0
# t_cell_data <- t_cell_data[,x]
# 
# t_cell_cor <- cor(t_cell_data[,1:sum(x)], method = 'spearman')
# 
# t_cell_cor_logic <- t_cell_cor > .7 & t_cell_cor < 1
# valid_n <- sum(t_cell_cor_logic)
# t_cell_cor_1 <- data.frame(row_gene = rep(0, times = valid_n), col_gene = rep(0, times = valid_n))
# n <- 1
# 
# for (i in 1:nrow(t_cell_cor)) {
#   for(j in 1:ncol(t_cell_cor)) {
#     if (t_cell_cor_logic[i,j] & !is.na(t_cell_cor_logic[i,j])) {
#       t_cell_cor_1$row_gene[n] <- rownames(t_cell_cor_logic)[i]
#       t_cell_cor_1$col_gene[n] <- colnames(t_cell_cor_logic)[j]
#       n <- n + 1
#     }
#   }
# }
# 
# t_cell_genes <- unique(t_cell_cor_1$row_gene)
# 
# # NK (> .7)
# nk_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[5]) %>% select(cell_type_gene_list[[5]], 'Simple_Relabelling', 'Patient_ID', 'State')
# 
# # Check gene standard deviation is > 0
# nk_sd <- apply(nk_data[,1:length(cell_type_gene_list[[5]])], 2, sd)
# x <- nk_sd > 0
# nk_data <- nk_data[,x]
# 
# nk_cor <- cor(nk_data[,1:sum(x)], method = 'spearman')
# 
# nk_cor_logic <- nk_cor > .7 & nk_cor < 1
# valid_n <- sum(nk_cor_logic)
# nk_cor_1 <- data.frame(row_gene = rep(0, times = valid_n), col_gene = rep(0, times = valid_n))
# n <- 1
# 
# for (i in 1:nrow(nk_cor)) {
#   for(j in 1:ncol(nk_cor)) {
#     if (nk_cor_logic[i,j] & !is.na(nk_cor_logic[i,j])) {
#       nk_cor_1$row_gene[n] <- rownames(nk_cor_logic)[i]
#       nk_cor_1$col_gene[n] <- colnames(nk_cor_logic)[j]
#       n <- n + 1
#     }
#   }
# }
# 
# nk_genes <- unique(nk_cor_1$row_gene)
# 
# # export each cluster out with filtered genes
# supplemental_cols <- c('Simple_Relabelling','Patient_ID','State')
# write.csv(myeloid_data[,c(myeloid_genes,supplemental_cols)], file = 'bioinf590_project_myeloid_filt.csv')
# write.csv(fibroblast_data[,c(fibroblast_genes,supplemental_cols)], file = 'bioinf590_project_fibroblast_filt.csv')
# write.csv(epithelial_data[,c(epithelial_genes,supplemental_cols)], file = 'bioinf590_project_epithelial_filt.csv')
# write.csv(t_cell_data[,c(t_cell_genes,supplemental_cols)], file = 'bioinf590_project_t_cell_filt.csv')
# write.csv(nk_data[,c(nk_genes,supplemental_cols)], file = 'bioinf590_project_nk_filt.csv')
# 
# '------------------------------------------------------------------------------------------------------------------------'
# # normalizing with z-scores to account for batch correction
# x <- tissue_data %>% filter(Patient_ID == '181429')
# x_scale <- scale(x[,1:4024], center = T, scale = T)
# y <- tissue_data %>% filter(Patient_ID == "19_227")
# y_scale <- scale(x[,1:4024], center = T, scale = T)
# z <- rbind(x_scale[1:26, 1:2],y_scale[1:72, 1:2])
# 
# raw_z <- rbind(x[1:26, 1:2], y[1:72, 1:2])
# 
# cor(raw_z, method = 'spearman')
# plot(z[1], z[2])
# 
# # Check zero inflation in each gene for each cell type
# x <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[1]) %>% select('SPP1')
# table(x > 0)
# y <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[1]) %>% select(cell_type_gene_list$Myeloid[2])
# table(y > 0)
# z <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[1]) %>% select(cell_type_gene_list$Myeloid[1:2])
# table(z > 0)
# 
# plot(tissue_data$SPP1[which(tissue_data$Simple_Relabelling == 'Myeloid')], tissue_data$S100A9[which(tissue_data$Simple_Relabelling == 'Myeloid')])
# cor(tissue_data$SPP1[which(tissue_data$Simple_Relabelling == 'Myeloid')], tissue_data$S100A9[which(tissue_data$Simple_Relabelling == 'Myeloid')])

'------------------------'
# 11/8/19

# load in data
load("~/rao_lab/human_interactome/Final_Tissue_Merge_Batch_Fil.RData")

# All genes in object
all_genes <- rownames(Tissue_Merge)

# Find common cell types for patients (Simple_Relabelling)
# Don't use 181658
# Use Myeloid, Epithelial, Fibroblast, T Cell, NK
table(Tissue_Merge$Simple_Relabelling, Tissue_Merge$Patient_ID) > 0
wanted_cell_type <- c('Myeloid', 'Epithelial', 'Fibroblast', 'T Cell', 'NK')

# Subset wanted cell types into new object
Idents(Tissue_Merge) <- 'Simple_Relabelling'
tissue_merge_class <- subset(Tissue_Merge, idents = wanted_cell_type)

# Run DE with batch correction between groups for each cell type 
tissue_merge_class_DE <- FindAllMarkers(tissue_merge_class, test.use = 'LR', latent.vars = 'Patient_ID')

# Filter DE genes of RP and MT- genes
filt_genes <- unique(tissue_merge_class_DE$gene[grep('^RP|^MT-', tissue_merge_class_DE$gene)])
DE_filt <- tissue_merge_class_DE %>% group_by(cluster) %>% filter(!gene %in% filt_genes)

# Get genes where p_val_adj < .05
DE_filt <- DE_filt %>% group_by(cluster) %>% filter(p_val_adj < .05)
cell_type_gene_list <- lapply(wanted_cell_type, function(x) DE_filt$gene[which(DE_filt$cluster == x)])
names(cell_type_gene_list) <- wanted_cell_type

# Create data frame of DE genes and wanted cell types
tissue_data <- FetchData(Tissue_Merge, vars = c(unique(DE_filt$gene), 'Simple_Relabelling', 'Patient_ID', 'State'))
tissue_data <- tissue_data %>% filter(Simple_Relabelling %in% wanted_cell_type & Patient_ID != '181658')

# Myeloid
## subset myeloid cells and DE genes
myeloid_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[1]) %>% select(cell_type_gene_list[[1]], 'Simple_Relabelling', 'Patient_ID', 'State')
myeloid_n_genes <- length(cell_type_gene_list[[1]])

## check for zero variance among myeloid DE genes
myeloid_var <- nearZeroVar(myeloid_data[,1:myeloid_n_genes], saveMetrics = T)
table(myeloid_var$zeroVar)

## calculate correlation
myeloid_cor <- cor(myeloid_data[,1:myeloid_n_genes], method = 'spearman')
diag(myeloid_cor) <- 0

## filter out genes that have correlation > .9
gene_remove <- findCorrelation(myeloid_cor, cutoff = .9)

if (length(gene_remove) > 0) {
  myeloid_genes <- rownames(myeloid_cor)[-gene_remove]
} else {
  myeloid_genes <- rownames(myeloid_cor)
}


# Epithelial
## subset epithelial cells and DE genes
epithelial_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[2]) %>% select(cell_type_gene_list[[2]], 'Simple_Relabelling', 'Patient_ID', 'State')
epithelial_n_genes <- length(cell_type_gene_list[[2]])

## check for zero variance among epithelial DE genes
epithelial_var <- nearZeroVar(epithelial_data[,1:epithelial_n_genes], saveMetrics = T)
table(epithelial_var$zeroVar)

## calculate correlation
epithelial_cor <- cor(epithelial_data[,1:epithelial_n_genes], method = 'spearman')
diag(epithelial_cor) <- 0

## filter out genes that have correlation > .9
gene_remove <- findCorrelation(epithelial_cor, cutoff = .9)

if (length(gene_remove) > 0) {
  epithelial_genes <- rownames(epithelial_cor)[-gene_remove]
} else {
  epithelial_genes <- rownames(epithelial_cor)
}

# Fibroblast
## subset fibroblast cells and DE genes
fibroblast_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[3]) %>% select(cell_type_gene_list[[3]], 'Simple_Relabelling', 'Patient_ID', 'State')
fibroblast_n_genes <- length(cell_type_gene_list[[3]])

## check for zero variance among fibroblast DE genes
fibroblast_var <- nearZeroVar(fibroblast_data[,1:fibroblast_n_genes], saveMetrics = T)
table(fibroblast_var$zeroVar)
fibroblast_data <- fibroblast_data %>% select(rownames(fibroblast_var)[which(!fibroblast_var$zeroVar)], 'Simple_Relabelling', 'Patient_ID', 'State')
fibroblast_n_genes <- ncol(fibroblast_data)-3

## calculate correlation
fibroblast_cor <- cor(fibroblast_data[,1:fibroblast_n_genes], method = 'spearman')
diag(fibroblast_cor) <- 0

## filter out genes that have correlation > .9
gene_remove <- findCorrelation(fibroblast_cor, cutoff = .9)

if (length(gene_remove) > 0) {
  fibroblast_genes <- rownames(fibroblast_cor)[-gene_remove]
} else {
  fibroblast_genes <- rownames(fibroblast_cor)
}

# T Cell
## subset T cell cells and DE genes
t_cell_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[3]) %>% select(cell_type_gene_list[[4]], 'Simple_Relabelling', 'Patient_ID', 'State')
t_cell_n_genes <- length(cell_type_gene_list[[4]])

## check for zero variance among t_cell DE genes
t_cell_var <- nearZeroVar(t_cell_data[,1:t_cell_n_genes], saveMetrics = T)
table(t_cell_var$zeroVar)
t_cell_data <- t_cell_data %>% select(rownames(t_cell_var)[which(!t_cell_var$zeroVar)], 'Simple_Relabelling', 'Patient_ID', 'State')
t_cell_n_genes <- ncol(t_cell_data)-3

## calculate correlation
t_cell_cor <- cor(t_cell_data[,1:t_cell_n_genes], method = 'spearman')
diag(t_cell_cor) <- 0

## filter out genes that have correlation > .9
gene_remove <- findCorrelation(t_cell_cor, cutoff = .9)

if (length(gene_remove) > 0) {
  t_cell_genes <- rownames(t_cell_cor)[-gene_remove]
} else {
  t_cell_genes <- rownames(t_cell_cor)
}

# NK
## subset nk cells and DE genes
nk_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[5]) %>% select(cell_type_gene_list[[5]], 'Simple_Relabelling', 'Patient_ID', 'State')
nk_n_genes <- length(cell_type_gene_list[[5]])

## check for zero variance among nk DE genes
nk_var <- nearZeroVar(nk_data[,1:nk_n_genes], saveMetrics = T)
table(nk_var$zeroVar)
nk_data <- nk_data %>% select(rownames(nk_var)[which(!nk_var$zeroVar)], 'Simple_Relabelling', 'Patient_ID', 'State')
nk_n_genes <- ncol(nk_data)-3

## calculate correlation
nk_cor <- cor(nk_data[,1:nk_n_genes], method = 'spearman')
diag(nk_cor) <- 0

## filter out genes that have correlation > .9
gene_remove <- findCorrelation(nk_cor, cutoff = .9)

if (length(gene_remove) > 0) {
  nk_genes <- rownames(nk_cor)[-gene_remove]
} else {
  nk_genes <- rownames(nk_cor)
}

# export with filtered genes
supplemental_cols <- c('Simple_Relabelling','Patient_ID','State')
write.csv(myeloid_data[,c(myeloid_genes,supplemental_cols)], file = 'bioinf590_project_myeloid_filt.csv')
write.csv(fibroblast_data[,c(fibroblast_genes,supplemental_cols)], file = 'bioinf590_project_fibroblast_filt.csv')
write.csv(epithelial_data[,c(epithelial_genes,supplemental_cols)], file = 'bioinf590_project_epithelial_filt.csv')
write.csv(t_cell_data[,c(t_cell_genes,supplemental_cols)], file = 'bioinf590_project_t_cell_filt.csv')
write.csv(nk_data[,c(nk_genes,supplemental_cols)], file = 'bioinf590_project_nk_filt.csv')


