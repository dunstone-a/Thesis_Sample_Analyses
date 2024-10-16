# Natri_2024_Lung: Analysis using fewer cell types
# Amelia Dunstone
# 2024-10-02

# Downsample the Natri (2024) lung data to have 2-3 cell types for each of the 
# 4 lineages. There are 10 states in total instead of the 38 that are used in the 
# paper. 


library(here)
library(dplyr)

# Require newest version for data to load correctly
# devtools::install_git("https://github.com/dunstone-a/QTLExperiment.git")
# BiocManager::install("multistateQTL")
library(QTLExperiment)
library(multistateQTL)

# Working directory is "Thesis_Sample_Analyses"

# Load LIMIX data -------------------------------------------------------------

input_path <- "Natri_2024_Lung/data/limix/"

state <- read.table(file.path(input_path, "filenames_limix"))$V1
state <- gsub("_qtl_results_all.txt", "", state)

input <- data.frame(
    state=state, 
    path=paste0(input_path, state, "_qtl_results_all.txt"))
input$broad_type <- sapply(strsplit(input$state, "_"), "[", 1)

# Most abundant cell types (most number of people with at least 5 cell types)

# Immune: MDM, Inflammatory monocyte, moDC,
# 
# Mesenchymal: Advential FB, Alveolar FB, SMC?
#   
# Epithelial: Alveolar type 2, Ciliated, Sec - SCGB1A1+/SCGB3A2,
# 
# Endothelial: Arteriole, Venule, aCap?

# 10 states, 11 minutes to load
selected_states <- c("immune_Monocyte-derivedmacrophage", 
    "immune_Inflammatorymonocyte", 
    "immune_moDC", 
    "mesenchymal_AdventitialFB",
    "mesenchymal_AlveolarFB",
    "epithelial_AT2",
    "epithelial_Ciliated",
    "epithelial_Secretory-SCGB1A1+SCGB3A2+",
    "endothelial_arteriole",
    "endothelial_venule")

input_reduced <- input[input$state %in% selected_states,]

start <- Sys.time()
# 11.86 minutes
qtle <- sumstats2qtle(
    input_reduced, 
    feature_id="feature_id",
    variant_id="snp_id", 
    betas="beta", 
    errors="beta_se",
    pvalues="p_value", 
    verbose=TRUE)
qtle
Sys.time() - start

betas(qtle) <- as.matrix(betas(qtle))
errors(qtle) <- as.matrix(errors(qtle))
pvalues(qtle) <- as.matrix(pvalues(qtle))

# Save LIMIX data --------------------------------------------------------------

# pretty quick, a few minutes
saveRDS(qtle, file = here("Natri_2024_Lung/data/qtle_limix_10_states.rds"))



# Load mashr data --------------------------------------------------------------

# No longer using this because
# a) For some reason, the mashr data had duplicate values for some LFSRs in some 
#    rows. This meant that operations like pvalues(qtle) can't be performed unless
#    I take the average of these elements. I modified sumstats2qtle but that made 
#    the function slower in all cases even when there are no duplicate values. 
# b) I am running mashr myself now in the script mashr_eQTL_analysis.R. 

input_path <- "Natri_2024_Lung/data/mashr/"
state <- read.table(paste0(input_path, "filenames_mashr"))$V1
state <- gsub(".tsv", "", state)

input <- data.frame(
    state=state,
    path=paste0(input_path, state, ".tsv"))
input$broad_type <- sapply(strsplit(input$state, "_"), "[", 1)

# Subset to the 10 states selected earlier
input_reduced <- input[input$state %in% selected_states,]


# mashr
# I dont know why but this only works when I split it up into smaller pieces! 
# I think it takes up too much RAM? I am using 64GB RAM 
# Times: 9 mins, 15 minutes, 12 minute, 9 minutes
# Total time is around 50 minutes. 
start <- Sys.time()
qtle_endo <- sumstats2qtle(
    input_reduced[1:2, ],
    feature_id="feature_id",
    variant_id="snp_rsid",
    betas="posterior_means",
    errors="sds",
    pvalues="lfsr",
    verbose=TRUE)
qtle_epi <- sumstats2qtle(
    input_reduced[3:5, ],
    feature_id="feature_id",
    variant_id="snp_rsid",
    betas="posterior_means",
    errors="sds",
    pvalues="lfsr",
    verbose=TRUE)
qtle_imm <- sumstats2qtle(
    input_reduced[6:8, ],
    feature_id="feature_id",
    variant_id="snp_rsid",
    betas="posterior_means",
    errors="sds",
    pvalues="lfsr",
    verbose=TRUE)
qtle_mes <- sumstats2qtle(
    input_reduced[9:10, ],
    feature_id="feature_id",
    variant_id="snp_rsid",
    betas="posterior_means",
    errors="sds",
    pvalues="lfsr",
    verbose=TRUE)
Sys.time() - start

# These all have the exact same rownames
all(
    rownames(qtle_endo) == rownames(qtle_epi),
    rownames(qtle_endo) == rownames(qtle_imm),
    rownames(qtle_endo) == rownames(qtle_mes))


# Choose 1 million rows to keep, about 10% of the data
ds <- sample(nrow(qtle_endo), size = 1000000)

qtle_mashr_ds <- cbind(
    qtle_endo[ds, ],
    qtle_epi[ds, ],
    qtle_imm[ds, ], 
    qtle_mes[ds, ]
)

# Changing the name of the test statistic assay
lfsrs(qtle_mashr_ds) <- pvalues(qtle_mashr_ds)
pvalues(qtle_mashr_ds) <- NULL

# Save mashr data --------------------------------------------------------------

# Save downsampled data
saveRDS(qtle_mashr_ds, file = here("Natri_2024_Lung/data/qtle_mashr_ds_10_states.rds"))

# Combine together the QTLE objects
# It was quicker to column bind this data than to select a subset of rows
start <- Sys.time()
qtle_mashr <- cbind(qtle_endo, qtle_epi, qtle_imm, qtle_mes)
Sys.time() - start

# Changing the name of the test statistic assay
lfsrs(qtle_mashr) <- pvalues(qtle_mashr)
pvalues(qtle_mashr) <- NULL


# Save full data
saveRDS(qtle_mashr, file = here("Natri_2024_Lung/data/qtle_mashr_10_states.rds"))

# mashr files have 6 columns and 15 million rows. 
# limix files have 12 million rows



