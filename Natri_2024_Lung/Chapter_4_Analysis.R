# Natri_2024_Lung: Analysis using fewer cell types
# Amelia Dunstone
# 2024-08-07

# Two Options:
# 8-10 cell types across 4 lineages
# 3 cell types across 3 lineages

library(here)

# BiocManager::install("QTLExperiment")
# BiocManager::install("multistateQTL")
library(QTLExperiment)
library(multistateQTL)

setwd(here("Natri_2024_Lung"))

# Load data --------------------------------------------------------------------

input_path <- "data/limix/"

state <- read.table("data/limix/filenames_limix")$V1
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

# 12 states
selected_states <- c("immune_Monocyte-derivedmacrophage", 
    "immune_Inflammatorymonocyte", 
    "immune_moDC", 
    "mesenchymal_AdventitialFB",
    "mesenchymal_AlveolarFB",
    "mesenchymal_SMC",
    "epithelial_AT2",
    "epithelial_Ciliated",
    "epithelial_Secretory-SCGB1A1+SCGB3A2+",
    "endothelial_arteriole",
    "endothelial_venule",
    "endothelial_aCap")

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

# 11 mins to load 
start <- Sys.time()
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

# pretty quick, a few minutes
saveRDS(qtle, file = "data/qtle_10_states.rds")

# QTLExperiment Vignette -------------------------------------------------------

## 3. Basic object manipulation 

dim(qtle)
# [1] 13941180       10
colnames(qtle)
head(rowData(qtle))

state_id(qtle) <- sapply(strsplit(state_id(qtle), "_"), "[", 2)
# Avoiding sapply:
state_id(qtle) <- vapply(strsplit(state_id(qtle), "_"), "[", character(1), 2)

# 4. Working with assays

# Getters and setters
betas(qtle)[1:5,1:5]

# 5. Working with critical meta data 

state_id(qtle)

feature_id(qtle)[1:3]

head(betas(qtle))

# MultistateQTL Vignette -------------------------------------------------------

## 3. Dealing with missing data ------------------------------------------------

devtools::load_all("~/Scratch/multistateQTL")
library(multistateQTL)

table(apply(is.na(pvalues(qtle)), MARGIN = 1, FUN = sum))
# 0       1       2       3       4       5       6       7       8 
# 8606050 3616167  465001  270701  228718  217378  196454  159118  114777 
# 9 
# 66816 

# For now doing a strict subset. 
# Can only have a maximum of 1 NA value for each test.
qtle_na <- getComplete(qtle, n=0.2, verbose=TRUE)

# Then for the remaining QTL, we can fill in the missing values using the 
# following scheme

# Replacing NAs doesn't work if using n = 0.1, because there are NO NAs
qtle_na <- replaceNAs(qtle_na)

head(betas(qtle_na))

## 4. Calling significance ------------------------------------------------------

# or 0.001? I feel like there needs to be some info about good thresholds to use. 
# Q for Davis?
qtle_na <- callSignificance(qtle_na, assay="pvalues", thresh=0.001)

message("Median number of significant tests per state: ", 
        median(colData(qtle_na)$nSignificant))
# Median number of significant tests per state: 18999.5

## 5. Plotting global patterns of sharing --------------------------------------

qtle_sig <- getSignificant(qtle_na)
qtle_top <- getTopHits(qtle_na, assay="pvalues", mode="state")
qtle_top <- runPairwiseSharing(qtle_top)

dim(qtle_sig)
dim(qtle_top)

p1 <- plotPairwiseSharing(qtle_top, annotate_cols=c("nSignificant", "broad_type"))
p1

# There is more sharing between the immune cell types than between the other 
# broad cell types, based on the number of top eQTLs which are shared. 


plotUpSet(qtle_top, annotateColsBy=c("nSignificant", "broad_type"), set_order = state_order)


state_order <- vapply(strsplit(selected_states, "_"), "[", character(1), 2)

## 6. Characterising multistateQTL patterns

qtle_top <- runTestMetrics(qtle_top)

plotCompareStates(qtle_top, x="venule", y="arteriole")

plotCompareStates(qtle_top, x="Inflammatorymonocyte", y="moDC")

# remove not sigs?

table(rowData(qtle_top)$qtl_type)
## 
##     global_diverging        global_shared multistate_diverging 
##                  253                  124                   63 
##    multistate_shared 
##                   13
hist(rowData(qtle_top)$nSignificant)


qtle_top_ms <- subset(qtle_top, qtl_type_simple == "multistate")

# qtle_top_ms <- subset(qtle_top_ms, )
plotQTLClusters(
  qtle_top_ms, 
  annotate_states="broad_type",
  annotate_tests=c("qtl_type"))

qtle_top_unique <- subset(qtle_top, qtl_type_simple == "unique")

# qtle_top_ms <- subset(qtle_top_ms, )
plotQTLClusters(
  qtle_top_unique, 
  annotate_states="broad_type",
  annotate_tests=c("qtl_type"))

