#!/usr/bin/env Rscript

# Network-based clustering of mutational and covariate data
library(clustNet)

# load data
mut_cov_data <- read.table("../data/categorical-mutationCovariate-matrix.txt")

mut_cov_data <- mut_cov_data[,-c(39:60)] # remove age, sex, cancer type

# clustering
startseed <- 100
nIterations <- 50
chi <- chi*0.1

bestCluster <- clustNet:::BBMMclusterEM(binaryMatrix = mut_cov_data,
                            chi = chi, k_clust = k,
                            startseed = startseed,
                            nIterations = nIterations)

saveRDS(bestCluster, paste0("euler_results/aic_k", k, "_chi", chi, ".rds"))


