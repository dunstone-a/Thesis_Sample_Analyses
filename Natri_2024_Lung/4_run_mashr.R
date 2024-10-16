# 'eQTL analysis outline' pipeline
# Amelia Dunstone
# 2024-10-01

# Following the pipeline to analyse eQTL data with mashr:
# https://stephenslab.github.io/mashr/articles/eQTL_outline.html 
# This is a modified version of the pipeline that is tailored for eQTL data. 
# As there are millions of tests in eQTL data, smaller subsets of the data are 
# used for particular steps.
# I have used QTLExperiment functionality to obtain the "strong" subset. 

# - 'strong' subset: Contains a subset of association tests where the results are
# more significant. This might be "top" eQTL in each gene based on univariate test 
# results e.g. "ashr". Urbut used 16,000 genes as their top subset (top eQTLs for each gene)
# - 'random' subset: A random selection of tests that should be unbiased. Should 
# include null and non-null tests. Urbut used 20,000 random tests. 

# 0.1 Load data -----------------------------------------------------------------

library(here)
library(QTLExperiment)
library(multistateQTL)

# Load data
# qtle <- readRDS(file = here("Natri_2024_Lung/data/qtle_limix_10_states.rds"))
root <- "/mnt/beegfs/mccarthy/scratch/general/adunstone/Thesis_Sample_Analyses"
qtle <- readRDS(file = file.path(root, "Natri_2024_Lung/data/qtle_limix_10_states.rds"))

message("Loaded LIMIX data as a QTLE")


# 0.2. Create random subset ----------------------------------------------------

# Identify a random subset of 5000 tests
random.subset <- sample(nrow(qtle), 5000)
length(random.subset)
# [1] 5000

# 0.3. Create strong subset ----------------------------------------------------

# Strong subject is the top eQTL for each gene

# Remove and impute NA values
qtle_na <- getComplete(qtle, n=0.5, verbose=TRUE)
qtle_na <- replaceNAs(qtle_na)

# Take the top test for each gene for the 'strong' subset
# 9118 genes
qtle_strong <- getTopHits(qtle_na, assay="pvalues", mode="global")

# Convert to just indices
strong.subset <- which(rownames(qtle) %in% rownames(qtle_strong))
length(strong.subset)
# [1] 9118


# 0.4. Convert to mashr object -------------------------------------------------

library(ashr)
library(mashr)
set.seed(1)

# Convert to mash object
mashed <- mash_set_data(Bhat = betas(qtle), Shat = errors(qtle))

# 1. Learn correlation structure -----------------------------------------------'

# Learn correlation structure among null tests using random test.

data.temp <- mash_set_data(mashed$Bhat[random.subset,], mashed$Shat[random.subset,])
Vhat <- estimate_null_correlation_simple(data.temp)
rm(data.temp)

data.random <- mash_set_data(mashed$Bhat[random.subset,], mashed$Shat[random.subset,], V=Vhat)
data.strong <- mash_set_data(mashed$Bhat[strong.subset,], mashed$Shat[strong.subset,], V=Vhat)

message("Created strong and random subsets")


# 2. Learn data-driven covariance matrices -------------------------------------

# Learn data-driven covariance matrices using strong tests. 

U.pca <- cov_pca(data.strong, 5)
message("Created U.pca")
# I get a warning:
# OpenBLAS Warning: Detect OpenMP Loop and this application may hang. Please 
# rebuild the library with USE_OPENMP=1 option. 

# NOTE:
# This piece of code takes a very long time to run. 
# Submitted as a bash script with 48h, see run_mashr.sh
U.ed <- suppressWarnings(cov_ed(data.strong, U.pca))
message("Created U.ed")


# 3. Fit the mashr model -------------------------------------------------------

# Fit the mashr model to the random tests, to learn the mixture weights on all 
# the different covariance matrices and scaling coefficients.

U.c <- cov_canonical(data.random)

m <- mash(data.random, Ulist = c(U.ed, U.c), outputlevel = 1)
#  - Computing likelihood matrices
message("Created m")

saveRDS(m, file = file.path(root, "Natri_2024_Lung/data/mashr_10_states.rds"))

message("Saved m")


# The outputlevel=1 option means that it will not compute posterior summaries for these tests (which saves time).



# 4. Compute posterior summaries -----------------------------------------------
# Compute posterior summaries on the strong tests, using the model fit from step 2. 
# (At this stage you could actually compute posterior summaries for any sets of tests 
# you like. For example you could read in all your tests in small batches and compute 
# posterior summaries in batches. But for illustration we will just do it on the strong tests.)


# (In mash the parameter g is used to denote the mixture model which we learned above.)

m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
#  - Computing 1428 x 241 likelihood matrix.
#  - Likelihood calculations took 0.07 seconds.
#  - Computing posterior matrices.
#  - Computation allocated took 0.03 seconds.
head(get_lfsr(m2))
#               condition_1  condition_2  condition_3  condition_4  condition_5
# effect_13096 9.815945e-06 5.056808e-01 4.229107e-01 3.944224e-01 6.055467e-01
# effect_29826 6.571537e-05 6.637417e-01 5.837333e-01 6.358124e-01 5.768253e-01
# effect_14042 6.994353e-02 6.495479e-03 2.483348e-03 5.562270e-02 6.836385e-06
# effect_12524 1.119195e-01 4.107543e-01 2.985565e-02 2.579205e-05 1.001824e-01
# effect_15456 4.913414e-05 4.380260e-01 2.733414e-01 5.166882e-01 3.610422e-01
# effect_35844 2.623221e-09 4.570036e-09 1.864892e-07 1.013875e-09 4.094924e-11


# Test:
tmp = simple_sims(10000,5,1) # simulates data on 40k tests
qtle <- mash2qtle(tmp, sep="\\|")

