# This script creates a figure for AIC analysis, used to determine the optimal number of clusters.

# Load required libraries
library(ggplot2)
library(reshape2)
library(clustNet)
library(extrafont)

# Clear the workspace
rm(list=ls())

# Helper Function to Load AIC Matrix
load_aic_matrix <- function(cluster_range, chi_range, file_path) {
  aic_matrix <- matrix(NA, nrow = length(cluster_range), ncol = length(chi_range))
  for (k_idx in seq_along(cluster_range)) {
    k <- cluster_range[k_idx]
    for (i_idx in seq_along(chi_range)) {
      i <- chi_range[i_idx]
      file <- paste0(file_path, "aic_k", k, "_chi", i, ".rds")
      if (file.exists(file)) {
        aic_value <- readRDS(file)$testAIC
        if (length(levels(as.factor(readRDS(file)$newclustermembership))) == k) {
          aic_matrix[k_idx, i_idx] <- aic_value
        }
      }
    }
  }
  return(aic_matrix)
}

# Helper Function to Plot AIC Heatmap
plot_aic_heatmap <- function(aics, minK, maxK, chiVec, AICrange) {
  meltdivergy <- melt(aics)
  ggheatmap <- ggplot(data = meltdivergy, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    xlab(expression(chi)) +
    ylab("k") +
    scale_fill_gradient2(high = rgb(0.98, 0.98, 1), low = "#117777",
                         mid = "#88BBBB", space = "Lab", na.value = "grey75",
                         midpoint = AICrange / 2, limit = c(0, AICrange), name = "AIC\nchange\n") +
    scale_y_continuous(breaks = c(minK:maxK)) +
    theme_minimal() +
    theme(axis.title.x = element_text(vjust = -1),
          axis.title.y = element_text(angle = 0, hjust = -0.5, vjust = 0.505),
          axis.text.x = element_text(angle = 0, vjust = 0.5, size = 20, hjust = 0.6),
          axis.text.y = element_text(angle = 0, vjust = 0.5, size = 20, hjust = 1),
          legend.text = element_text(size = 20),
          axis.title = element_text(size = 30),
          legend.title = element_text(size = 24)) +
    theme(legend.key.size = unit(2, "line")) +
    theme(plot.margin = unit(c(-0.3, -0.3, 0.4, 0.4), "cm"))
  return(ggheatmap)
}

# Main Script

# Initialize Parameters
cluster_range <- 5:18
chi_range <- 0:40 * 0.1
file_path <- "../euler_AIC/euler_results/"

# Load AIC Matrix
aic_matrix <- load_aic_matrix(cluster_range, chi_range, file_path)

# Prepare and Plot AIC Heatmap
minK <- 5
maxK <- 18
AICrange <- 40
minaics <- apply(aic_matrix, 2, min, na.rm = TRUE)
aics <- t(aic_matrix) - minaics
topaics <- AICrange
aics[aics > topaics] <- topaics
rownames(aics) <- chi_range
colnames(aics) <- c(minK:maxK)

ggheatmap <- plot_aic_heatmap(aics, minK, maxK, chi_range, AICrange)

# Save Plot to PDF
loadfonts()
pdf("~/Desktop/aic_analysis2.pdf", height = 8, width = 16, family = "Arial", paper = "special", onefile = FALSE)
ggheatmap
dev.off()


