# scRNA Covariance Patient Classification
# 11/4/19

library(dplyr)
library(tibble)
library(Seurat)
library(corrplot)

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

# Get genes where p_val_adj < .05
DE_filt <- tissue_merge_class_DE %>% group_by(cluster) %>% filter(p_val_adj < .05)
cell_type_gene_list <- lapply(wanted_cell_type, function(x) DE_filt$gene[which(DE_filt$cluster == x)])
names(cell_type_gene_list) <- wanted_cell_type

# Create data frame of DE genes and wanted cell types
tissue_data <- FetchData(Tissue_Merge, vars = c(unique(DE_filt$gene), 'Simple_Relabelling', 'Patient_ID', 'State'))
tissue_data <- tissue_data %>% filter(Simple_Relabelling %in% wanted_cell_type & Patient_ID != '181658')

# Myeloid (> .8)
myeloid_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[1]) %>% select(cell_type_gene_list[[1]], 'Simple_Relabelling', 'Patient_ID', 'State')
myeloid_cor <- cor(myeloid_data[,1:length(cell_type_gene_list[[1]])], method = 'spearman')
table(abs(myeloid_cor) > .7 & abs(myeloid_cor) != 1)
corrplot(abs(myeloid_cor[1:30, 1:30]) > .7 & abs(myeloid_cor[1:30, 1:30]) < 1)

# Epithelial (> .7)
epithelial_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[2]) %>% select(cell_type_gene_list[[2]], 'Simple_Relabelling', 'Patient_ID', 'State')
epithelial_cor <- cor(epithelial_data[,1:length(cell_type_gene_list[[2]])], method = 'spearman')
table(abs(epithelial_cor) > .7 & abs(epithelial_cor) != 1)
corrplot(abs(epithelial_cor[1:30, 1:30]) > .7 & abs(epithelial_cor[1:30, 1:30]) < 1)

# Fibroblast (> .9)
fibroblast_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[3]) %>% select(cell_type_gene_list[[3]], 'Simple_Relabelling', 'Patient_ID', 'State')
fibroblast_cor <- cor(fibroblast_data[,1:length(cell_type_gene_list[[3]])], method = 'spearman')
table(abs(fibroblast_cor) > .7 & abs(fibroblast_cor) != 1)
corrplot(abs(fibroblast_cor[1:30, 1:30]) > .7 & abs(fibroblast_cor[1:30, 1:30]) < 1)

# T Cell (> .8)
t_cell_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[4]) %>% select(cell_type_gene_list[[4]], 'Simple_Relabelling', 'Patient_ID', 'State')
t_cell_cor <- cor(t_cell_data[,1:length(cell_type_gene_list[[4]])], method = 'spearman')
table(abs(t_cell_cor) > .7 & abs(t_cell_cor) != 1)
corrplot(abs(t_cell_cor[1:30, 1:30]) > .7 & abs(t_cell_cor[1:30, 1:30]) < 1)

# NK (> .7)
nk_data <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[5]) %>% select(cell_type_gene_list[[5]], 'Simple_Relabelling', 'Patient_ID', 'State')
nk_cor <- cor(nk_data[,1:length(cell_type_gene_list[[5]])], method = 'spearman')
table(abs(nk_cor) > .7 & abs(nk_cor) != 1)
corrplot(abs(nk_cor[1:30, 1:30]) > .7 & abs(nk_cor[1:30, 1:30]) < 1)

'------------------------------------------------------------------------------------------------------------------------'
# normalizing with z-scores to account for batch correction
x <- tissue_data %>% filter(Patient_ID == '181429')
x_scale <- scale(x[,1:4024], center = T, scale = T)
y <- tissue_data %>% filter(Patient_ID == "19_227")
y_scale <- scale(x[,1:4024], center = T, scale = T)
z <- rbind(x_scale[1:26, 1:2],y_scale[1:72, 1:2])

raw_z <- rbind(x[1:26, 1:2], y[1:72, 1:2])

cor(raw_z, method = 'spearman')
plot(z[1], z[2])

# Check zero inflation in each gene for each cell type
x <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[1]) %>% select('SPP1')
table(x > 0)
y <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[1]) %>% select(cell_type_gene_list$Myeloid[2])
table(y > 0)
z <- tissue_data %>% filter(Simple_Relabelling == wanted_cell_type[1]) %>% select(cell_type_gene_list$Myeloid[1:2])
table(z > 0)

plot(tissue_data$SPP1[which(tissue_data$Simple_Relabelling == 'Myeloid')], tissue_data$S100A9[which(tissue_data$Simple_Relabelling == 'Myeloid')])
cor(tissue_data$SPP1[which(tissue_data$Simple_Relabelling == 'Myeloid')], tissue_data$S100A9[which(tissue_data$Simple_Relabelling == 'Myeloid')])

