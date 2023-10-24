# Code for the summary figure

# Load required libraries
library(ggplot2)
library(ggpubr)
library(reshape2)

# Clear the workspace
rm(list=ls())

p1 <- readRDS("../figures/density_plot_cluster_ct_ct_clean.rds")
p2 <- readRDS("../figures/density_plot_cluster_cluster.rds")
p3 <- readRDS("../figures/barplot.rds")
p4 <- readRDS("../figures/km_plot.rds")

p_mc <- readRDS("../figures/mut_cyto_barplot.rds")

cluster_legend <- readRDS("../figures/bar_legend.rds")
cluster_legend_ct <- readRDS("../figures/bar_legend_ct.rds")

ggarrange(ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2,widths = c(1,1), heights = c(1,1), labels = c("a","b","c","d")), 
          ggarrange(cluster_legend_ct, cluster_legend), ncol = 1, nrow = 2, heights = c(1,0.08))

ggarrange(ggarrange(cluster_legend_ct, cluster_legend), ggarrange(ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2,widths = c(1,1), heights = c(1,1), labels = c("a","b","c","d")), 
                                                                  p_mc, ncol = 1, nrow = 2, heights = c(2,1.15), labels = c("", "e")), ncol = 1, nrow = 2, heights = c(0.08,1))

pdf("~/Desktop/density_panel2.pdf", width = 9, height = 7, onefile=FALSE)
ggarrange(ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2,widths = c(1,1), heights = c(1,1), labels = c("a","b","c","d")), 
          ggarrange(cluster_legend_ct, cluster_legend), ncol = 1, nrow = 2, heights = c(1,0.08))
dev.off()

pdf("~/Desktop/density_panel3.pdf", width = 9, height = 9, onefile=FALSE)
ggarrange(ggarrange(cluster_legend_ct, cluster_legend), ggarrange(ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2,widths = c(1,1), heights = c(1,1), labels = c("a","b","c","d")),
          p_mc, ncol = 1, nrow = 2, heights = c(2,1.15), labels = c("", "e")), ncol = 1, nrow = 2, heights = c(0.08,1))
dev.off()

