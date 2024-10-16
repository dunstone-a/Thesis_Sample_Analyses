# Subset the seurat object from Natri_2024_Lung
# Amelia Dunstone
# 2024-10-09

# Step 1: Subset Seurat object to 10 cell types.

# I am using 256GB RAM for this script, the Seurat object is 64GB of memory. 


# install.packages("SeuratObject")
# install.packages("Seurat")

library(here)
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(dplyr)

# 0. Load data -----------------------------------------------------------------

# Script to determine which cell types are the most abundant in the Natri (2024) 
# data. 


# Removed .gz extension using bash script `2_extract_data.sh`

# Root directory is Thesis_Sample_Analyses
root <- "/mnt/beegfs/mccarthy/scratch/general/adunstone/Thesis_Sample_Analyses"

# Load RDS file
lung <- readRDS(file = file.path(root, "Natri_data/seurat/Seurat.rds"))

print(class(lung))


# 1. Subset to 10 cell types ---------------------------------------------------


# table(donor, cell_type)
# rowSums
# See what the range is for how many cell types there are for each donor. 
print(table(lung@meta.data$Short_Sample_Name))
print(length(table(lung@meta.data$Short_Sample_Name)))

print(table(lung@meta.data$manual_annotation_1))
print(length(table(lung@meta.data$manual_annotation_1)))

print(rowSums(table(lung@meta.data$Short_Sample_Name, lung@meta.data$manual_annotation_1) >= 5))


selected_states <- c("Monocyte-derived macrophage", "Inflammatory monocyte", 
    "moDC", "Adventitial FB", "Alveolar FB", "AT2", "Ciliated",
    "Secretory - SCGB1A1+/SCGB3A2+", "Arteriole", "Venule")



lung_10_states <- subset(x = lung, subset = manual_annotation_1 %in% selected_states)

saveRDS(lung_10_states, file = file.path(root, "Natri_data/seurat/Seurat_10_states.rds"))

message("Saved object")
