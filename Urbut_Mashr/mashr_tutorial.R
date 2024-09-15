# mashr vignette
# Amelia Dunstone
# 2024-09-02

# Following tutorial in the mashr vignette. 
# Based on paper Urbut et al. 2018.

# There are four steps to a mashr analysis

# - Read in the data
# - Set up the covariance matrices to be used
# - Fit the model
# - Extract posterior summaries and other quantities of interest

# 'The crucial rule'
# - mashr must be performed on all tests you perform, or a random subset of tests. 
# It CANNOT be used for just the most "significant" tests. 
# This is because mashr tries to estimate the null signal in the data, similar to 
# a multiple testing correction. 

library(ashr)
library(mashr)
set.seed(1)

# It creates 500 effects of each type for a total of 2000 effects.
simdata = simple_sims(500,5,1)

## Step 1: Read in the data ----------------------------------------------------

data = mash_set_data(simdata$Bhat, simdata$Shat)

# "[If you have only access to Z scores, you can set Bhat to the Z scores, and 
# set Shat to be the matrix with all 1s]."


## Step 2: Set up the covariance matrices --------------------------------------

U.c = cov_canonical(data)  
print(names(U.c))

# Identity:
# Condition_i: e_i x e_i (matrix with one 1 on the i'th diagonal)
# equal effects: matrix of all 1s
# simple_het_1: diagonal of 1s, covariances of 0.25
# simple_het_2: diagonal of 1s, covariances of 0.5
# simple_het_3: diagonal of 1s, covariances of 0.75

## Step 3: fit the model -------------------------------------------------------

# grid parameter: NULL
m.c = mash(data, U.c)

## Step 4: Extract Posterior Summaries -----------------------------------------

# This step must be peformed using all the tests (or a large random subset)
# These are J x R matrices. (tests by conditions)

# Local false sign rate
head(get_lfsr(m.c))

# posterior mean 
head(get_pm(m.c))

# posteriore standard deviation
head(get_psd(m.c))

# Find significant tests, default threshold is 0.05.
# Retains tests that are significant in at least one condition. 
head(get_significant_results(m.c))
print(length(get_significant_results(m.c)))

# Get the significant results in a subset of conditions (just condition 1)
print(head(get_significant_results(m.c, conditions=1)))

## Sharing ---------------------------------------------------------------------

# sharing of significant signals among each pair of conditions
# Default value for sharing is: the same sign and within a factor of 0.5 of each other. 
print(get_pairwise_sharing(m.c)) 

# Use a factor of 0 to test only if the sign is the same 
print(get_pairwise_sharing(m.c, factor=0))

## # mashr vignette
# Amelia Dunstone
# 2024-09-02

# Following tutorial in the mashr vignette. 
# Based on paper Urbut et al. 2018.

# There are four steps to a mashr analysis

# - Read in the data
# - Set up the covariance matrices to be used
# - Fit the model
# - Extract posterior summaries and other quantities of interest

# 'The crucial rule'
# - mashr must be performed on all tests you perform, or a random subset of tests. 
# It CANNOT be used for just the most "significant" tests. 
# This is because mashr tries to estimate the null signal in the data, similar to 
# a multiple testing correction. 

library(ashr)
library(mashr)
set.seed(1)

# It creates 500 effects of each type for a total of 2000 effects.
simdata = simple_sims(500,5,1)

## Step 1: Read in the data ----------------------------------------------------

data = mash_set_data(simdata$Bhat, simdata$Shat)

# "[If you have only access to Z scores, you can set Bhat to the Z scores, and 
# set Shat to be the matrix with all 1s]."


## Step 2: Set up the covariance matrices --------------------------------------

U.c = cov_canonical(data)  
print(names(U.c))

# Identity:
# Condition_i: e_i x e_i (matrix with one 1 on the i'th diagonal)
# equal effects: matrix of all 1s
# simple_het_1: diagonal of 1s, covariances of 0.25
# simple_het_2: diagonal of 1s, covariances of 0.5
# simple_het_3: diagonal of 1s, covariances of 0.75

## Step 3: fit the model -------------------------------------------------------

# grid parameter: NULL
m.c = mash(data, U.c)

## Step 4: Extract Posterior Summaries -----------------------------------------

# This step must be peformed using all the tests (or a large random subset)
# These are J x R matrices. (tests by conditions)

# Local false sign rate
head(get_lfsr(m.c))

# posterior mean 
head(get_pm(m.c))

# posteriore standard deviation
head(get_psd(m.c))

# Find significant tests, default threshold is 0.05.
# Retains tests that are significant in at least one condition. 
head(get_significant_results(m.c))
print(length(get_significant_results(m.c)))

# Get the significant results in a subset of conditions (just condition 1)
print(head(get_significant_results(m.c, conditions=1)))

## Sharing ---------------------------------------------------------------------

# sharing of significant signals among each pair of conditions
# Default value for sharing is: the same sign and within a factor of 0.5 of each other. 
print(get_pairwise_sharing(m.c)) 

# Use a factor of 0 to test only if the sign is the same 
print(get_pairwise_sharing(m.c, factor=0))

## Measure of fit (Likelihood) -------------------------------------------------

print(get_loglik(m.c))

## Estimated mixture proportions -----------------------------------------------

print(get_estimated_pi(m.c))

# Here we can see most of the mass is on the null, identity, singletons_1 (which 
# corresponds to effects that are specific to condition 1) and equal_effects.
# This matches the way the data was generated. 

barplot(get_estimated_pi(m.c),las = 2)

## Metaplot --------------------------------------------------------------------

# This plot shows the effects in each condition for the most significant effect. 
# It can be seen that 
mash_plot_meta(m.c,get_significant_results(m.c)[1])
