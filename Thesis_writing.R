# Code to produce figures that are in the paper
# Amelia Dunstone
# 2024/10/10

# Linear model -----------------------------------------------------------------


set.seed(2)
library(ggbeeswarm)

n <- c(20, 15, 18)
mu <- c(4, 5.1, 6.2)
sigma <- 0.6
AA <- rnorm(n[1], mean = mu[1], sd = sigma)
Aa <- rnorm(n[2], mean = mu[2], sd = sigma)
aa <- rnorm(n[3], mean = mu[3], sd = sigma)
ls <- list(AA = AA, Aa = Aa, aa = aa)

boxplot(ls, xlab = "Genotype")


boxplot(ls, xlab = "s")










n <- 2
p <- c(0.48, 0.02)

E <- n*p
E
variance <- n*p*(1-p)
variance
