# Explore the seurat object from Natri_2024_Lung
# Amelia Dunstone
# 2024-10-09

# Step 2: UMAP plot of data and any other necessary scRNA-seq steps

# I am using 256GB RAM for this script, the Seurat object is 64GB of memory. 

library(here)
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(dplyr)

selected_states <- c("Monocyte-derived macrophage", "Inflammatory monocyte", 
    "moDC", "Adventitial FB", "Alveolar FB", "AT2", "Ciliated",
    "Secretory - SCGB1A1+/SCGB3A2+", "Arteriole", "Venule")

# Load data with 14 states --------------------------------------------------------------------

lung <- readRDS(here("Natri_2024_Lung/data/seurat/Seurat_10_states.rds"))

# This data has already been processed in the paper so there is no need to remove 
# any cells for quality control reasons, or correct batches etc. 
# Just running PCA and doing clustering and visualising. 

lung <- RunPCA(object = lung)
lung <- FindNeighbors(object = lung, dims = 1:30)
lung <- FindClusters(object = lung)
lung <- RunUMAP(object = lung, dims = 1:30)

num_groups <- nlevels(factor(lung$manual_annotation_1))

# UMAP coloured by cell types
DimPlot(object = lung,
    reduction = "umap",
    group.by = "manual_annotation_1",
    cols = RColorBrewer::brewer.pal(num_groups, "Set3"))

# UMAP coloured by lineage
DimPlot(object = lung,
        reduction = "umap",
        group.by = "lineage",
        cols = RColorBrewer::brewer.pal(num_groups, "Dark2"))

# Load full data (all states) ----------------------------------------------------------------


lung_all <- readRDS(here("Natri_2024_Lung/data/seurat/Seurat.rds"))

# This data has already been processed in the paper so there is no need to remove 
# any cells for quality control reasons, or correct batches etc. 
# Just running PCA and doing clustering and visualising. 

# PCA and run UMAP
lung_all <- RunPCA(object = lung_all)
lung_all <- FindNeighbors(object = lung_all, dims = 1:30)
lung_all <- FindClusters(object = lung_all)
lung_all <- RunUMAP(object = lung_all, dims = 1:30)

# Number of cell-types
num_groups <- nlevels(factor(lung_all$manual_annotation_1))

# UMAP coloured by all the 38 cell types
DimPlot(object = lung_all,
        reduction = "umap",
        group.by = "manual_annotation_1")


# Variable for cell-type but only the 14 that I am including
lung_all$in_analysis <- lung_all$manual_annotation_1
lung_all$in_analysis[!(lung_all$manual_annotation_1 %in% selected_states)] <- NA
lung_all$in_analysis <- factor(lung_all$in_analysis,
                               levels = selected_states)

# Selected cell-types with their lineage label
lung_all$lineage_in_analysis <- factor(lung_all$lineage,
                                       levels = c("Endothelial", "Epithelial", "Immune", "Mesenchymal"))
lung_all$lineage_in_analysis[!(lung_all$manual_annotation_1 %in% selected_states)] <- NA

celltype_colours <- setNames(nm = levels(lung_all$in_analysis),
    c(RColorBrewer::brewer.pal(4, "Set3")[c(1,1,)], "grey75"))

lineage_colours <- setNames(RColorBrewer::brewer.pal(4, "Dark2"),
    nm = c("Immune", "Mesenchymal", "Epithelial", "Endothelial"))

celltype_colours <- setNames(lineage_colours[c(
    rep("Immune", 3), 
    rep("Mesenchymal", 2), 
    rep("Epithelial", 3), 
    rep("Endothelial", 2))],
    nm = selected_states)

lung_all$short_names <- lung_all$in_analysis
levels(lung_all$short_names) <- c("MDM", "Inflammatory monocyte", "moDC",
    "Adventitial FB", "Alveolar FB", "AT2", "Ciliated", "Secretory",
    "Arteriole", "Venule")

short_colours <- setNames(celltype_colours[levels(lung_all$in_analysis)], 
                          nm = levels(lung_all$short_names))



p1 <- DimPlot(object = lung_all, 
        label = TRUE,
        reduction = "umap",
        group.by = "short_names",
        cols = short_colours,
        na.value = "grey75") + 
    ggtitle("UMAP plot of lung single-cell data")
p1

p2 <- DimPlot(object = lung_all, 
              label = TRUE,
              reduction = "umap",
              group.by = "lineage_in_analysis",
              cols = lineage_colours,
              na.value = "grey75") + 
    ggtitle("UMAP plot of lung single-cell data")
p2

png("Figures/scRNAseq.png", width = 2200, height = 1600, res = 300)
p1
dev.off()

png("Figures/scRNAseq_lineage.png", width = 2200, height = 1600, res = 300)
p2
dev.off()

