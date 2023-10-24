#!/usr/bin/env Rscript

# Network-based clustering of mutational data

library(netClust)

# load data
mutCovData <- read.table("../data/binary-mutationCovariate-matrix.txt")
stringEdgepmat <- as.matrix(read.table("../data/string-edgepmat.txt"))

# define covariates
Ngenes <- 201
Ncov <- length(colnames(mutCovData))-Ngenes
mutData <- mutCovData[,1:Ngenes]

# clustering
itLim <- 100
EMseeds <- 1
clusterResMut <- netClust(mutData, kclust = 22, nbg = 0,itLim = itLim, EMseeds = EMseeds, edgepmat = stringEdgepmat)
saveRDS(clusterResMut, "../results/BNMM_var_memberships.rds")

