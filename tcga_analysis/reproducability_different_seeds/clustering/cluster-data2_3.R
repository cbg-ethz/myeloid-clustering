#!/usr/bin/env Rscript

# Network-based clustering of mutational and covariate data

library(netClust)

# load data
mutCovData <- read.table("../data/binary-mutationCovariate-matrix.txt")
stringEdgepmat <- as.matrix(read.table("../data/string-edgepmat.txt"))

# define covariates
Nall <- length(colnames(mutCovData))
mutCovData <- cbind(mutCovData[1:201], mutCovData[204:Nall], mutCovData[202:203])
Ncov <- 2

# clustering
itLim <- 200
EMseeds <- 3
clusterResCond <- netClust(mutCovData, kclust = 22, nbg = Ncov,itLim = itLim, EMseeds = EMseeds, edgepmat = stringEdgepmat)
saveRDS(clusterResCond, paste0("clusterResCond_i",itLim, "_s", min(EMseeds), ".rds"))

