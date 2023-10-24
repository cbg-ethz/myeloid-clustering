# Network-based clustering of mutational and covariate data
library(clustNet)

# load mutationa and covariate data
mut_cov_data <- read.table("../data/binary-mutationCovariate-matrix.txt")[,-c(58:61)]

# clustering (might take a few minutes, alternatively just load the file in the next line)
cluster_results <- get_clusters(mut_cov_data, k_clust = 9, n_bg = 2, quick = FALSE) # n_cov defines the number of cluster-independent covariates (they need to be in the last cols of the matrix)

# or read if the clusters from file
cluster_results <- readRDS("../results/euler_memberships.rds")

# Load additional pacakges to visualize the networks
library(ggplot2)
library(ggraph)
library(igraph)
library(ggpubr)

# Visualize networks
plot_clusters(cluster_results, directed = FALSE, node_colours = c(rep("#92a8d1",38),rep("#064273",8),rep("#e56d78",6), rep("#f7cac9",5),  rep("#708090",2)))

# Load additional pacakges to create a 2d dimensionality reduction
library(car)
library(ks)
library(graphics)
library(stats)

# Plot a 2d dimensionality reduction
density_plot(cluster_results, var_selection = 1:52, colourys = c("#202020","#774411","#DDAA77","#ed2124","#114477","#CC99BB", "#88CCAA","#117744","#77AADD"))  # var_selection does focus on mutations + cytogenetics + blood values only
