# Download Natri_2024_Lung data from GEO
# Amelia Dunstone
# 2024-10-02

# Download Natri (2024) data from the link provided in the 'Data Availability' 
# section of the paper. 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227136

dir.create("Natri_2024_Lung/data")
dir.create("Natri_2024_Lung/data/seurat")
dir.create("Natri_2024_Lung/data/raw")
dir.create("Natri_2024_Lung/data/limix")
dir.create("Natri_2024_Lung/data/mashr")

# Download LIMIX and mashr data ------------------------------------------------

input_path <- data.frame(
    path = c("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE227136&format=file&file=GSE227136%5Flimix%5Fres%5Fendothelial%2Etar%2Egz",
             "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE227136&format=file&file=GSE227136%5Flimix%5Fres%5Fepithelial%2Etar%2Egz",
             "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE227136&format=file&file=GSE227136%5Flimix%5Fres%5Fimmune%2Etar%2Egz",
             "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE227136&format=file&file=GSE227136%5Flimix%5Fres%5Fmesenchymal%2Etar%2Egz",
             "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE227136&format=file&file=GSE227136%5Fmashr%5Fapplied%5Fendothelial%2Etar%2Egz",
             "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE227136&format=file&file=GSE227136%5Fmashr%5Fapplied%5Fepithelial%2Etar%2Egz",
             "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE227136&format=file&file=GSE227136%5Fmashr%5Fapplied%5Fimmune%2Etar%2Egz",
             "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE227136&format=file&file=GSE227136%5Fmashr%5Fapplied%5Fmesenchymal%2Etar%2Egz", 
             "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE227136&format=file"),
    type = c(rep(c("limix", "mashr"), each = 4), "raw"), 
    names = c("GSE227136_limix_res_endothelial.tar.gz",
              "GSE227136_limix_res_epithelial.tar.gz", 
              "GSE227136_limix_res_immune.tar.gz",
              "GSE227136_limix_res_mesenchymal.tar.gz", 
              "GSE227136_mashr_applied_endothelial.tar.gz",
              "GSE227136_mashr_applied_epithelial.tar.gz",
              "GSE227136_mashr_applied_immune.tar.gz", 
              "GSE227136_mashr_applied_mesenchymal.tar.gz", 
              "GSE227136_RAW.tar"))

input_path$dest <- file.path("Natri_2024_Lung/data", input_path$type, input_path$names)

# Increase timeout value
options(timeout=4000)

# Download data
download.file(input_path$path, input_path$dest, method="libcurl")

# Each of these tar files uncompresses to be multiple text files. 

# Best to use bash functions for extracting. See extract_data.sh

# Download Seurat object -------------------------------------------------------


# Increase timeout value
options(timeout=10000)

# Download data
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227136/suppl/GSE227136%5FILD%5Fall%5Fcelltypes%5FSeurat.rds.gz", 
  destfile = "Natri_2024_Lung/data/seurat/Seurat.rds.gz", 
  method ="libcurl")

# This file is 20GB so needed to increase timeout value. 




