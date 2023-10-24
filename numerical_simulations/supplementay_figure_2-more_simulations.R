
# Prepare session, load packages
rm(list=ls())
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(clustNet)
library(reshape2)
# library(extrafont)

#################################
## Varied Number of Covariates ##
#################################

## plot netClust versions for different numbers of covariates
res1 <- readRDS("results/results--k_clust-6--n_vars-20--n_bg-0--n_it-30--different--TRUE.rds")
# nice_plot(res1$correct_samples)
nbg_data <- melt(res1$correct_samples[,1:9])
nbg_data$nbg <- 0
res2 <- readRDS("results/results--k_clust-6--n_vars-20--n_bg-2--n_it-30--different--TRUE.rds")
# nice_plot(res2$correct_samples)
nbg_data_temp <- melt(res2$correct_samples[,1:9])
nbg_data_temp$nbg <- 2
nbg_data <- rbind(nbg_data,nbg_data_temp)
res3 <- readRDS("results/results--k_clust-6--n_vars-20--n_bg-4--n_it-30--different--TRUE.rds")
nbg_data_temp <- melt(res3$correct_samples[,1:9])
# nice_plot1 <- nice_plot(res3$correct_samples)
nbg_data_temp$nbg <- 4
nbg_data <- rbind(nbg_data,nbg_data_temp)
res4 <- readRDS("results/results--k_clust-6--n_vars-20--n_bg-6--n_it-30--different--TRUE.rds")
# nice_plot(res4$correct_samples)
nbg_data_temp <- melt(res4$correct_samples[,1:9])
nbg_data_temp$nbg <- 6
nbg_data <- rbind(nbg_data,nbg_data_temp)

colnames(nbg_data)[2] <- "Method"

# (publication figure 2)

# remove Mclust (since Gaussian mixture model not interesting for binary data)
nbg_data_reduced <- nbg_data[!nbg_data$Method=="Mclust (cov & var)"&!nbg_data$Method=="Mclust (var)",]
# create labels for publication
levels(nbg_data_reduced$Method) <- c("CANclust (cov. & var.)", "BNMM (cov. & var.)",
                             "BNMM (var.)", "K-means (cov. & var.)", "K-means (var.)",
                             "MC", "MC2", "BMM (cov. & var.)", "BMM (var.)")
# choose colors 
# color_list <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")
# cols <- color_list[c(2,9,4,1,10,5,3)]
# cols <- c("#E31A1C", "#114477", "#77AADD", "#771122","#DD7788", "#117744", "#88CCAA")
cols <- c("#E31A1C", "#117777", "#77CCCC", "#771122","#DD7788", "#777711", "#DDDD77")


plot_nbg <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted rand index (ARI)", x = "Number of covariates") +
  theme_minimal(); plot_nbg

# plot_nbg_nolegend <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
#   # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
#   geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
#   scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
#   labs(y = "Adjusted rand index (ARI)", x = "Number of covariates") +
#   theme_minimal()+
#   theme(legend.position = "none"); plot_nbg_nolegend

# # save
# png("~/Desktop/plot_nbg.png", width = 18, height = 7, units = 'cm', res = 300)
# plot_nbg
# dev.off()


###############################
## Varied Number of Clusters ##
###############################

## plot netClust versions for different numbers of covariates
res1 <- readRDS("results/results--k_clust-2--n_vars-20--n_bg-5--n_it-20--different--TRUE.rds")
# nice_plot(res1$correct_samples)
nbg_data <- melt(res1$correct_samples[,1:9])
nbg_data$nbg <- 2
res2 <- readRDS("results/results--k_clust-4--n_vars-20--n_bg-5--n_it-20--different--TRUE.rds")
# nice_plot(res2$correct_samples)
nbg_data_temp <- melt(res2$correct_samples[,1:9])
nbg_data_temp$nbg <- 4
nbg_data <- rbind(nbg_data,nbg_data_temp)
res3 <- readRDS("results/results--k_clust-6--n_vars-20--n_bg-5--n_it-20--different--TRUE.rds")
nbg_data_temp <- melt(res3$correct_samples[,1:9])
# nice_plot1 <- nice_plot(res3$correct_samples)
nbg_data_temp$nbg <- 6
nbg_data <- rbind(nbg_data,nbg_data_temp)
res4 <- readRDS("results/results--k_clust-8--n_vars-20--n_bg-5--n_it-20--different--TRUE.rds")
# nice_plot(res4$correct_samples)
nbg_data_temp <- melt(res4$correct_samples[,1:9])
nbg_data_temp$nbg <- 8
nbg_data <- rbind(nbg_data,nbg_data_temp)

colnames(nbg_data)[2] <- "Method"

# remove Mclust (since Gaussian mixture model not interesting for binary data)
nbg_data_reduced <- nbg_data[!nbg_data$Method=="Mclust (cov & var)"&!nbg_data$Method=="Mclust (var)",]
# create labels for publication
levels(nbg_data_reduced$Method) <- c("CANclust (cov. & var.)", "BNMM (cov. & var.)",
                                     "BNMM (var.)", "K-means (cov. & var.)", "K-means (var.)",
                                     "MC", "MC2", "BMM (cov. & var.)", "BMM (var.)")

# choose colors 
# color_list <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")
# cols <- color_list[c(2,9,4,1,10,5,3)]

plot_nclust <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted rand index (ARI)", x = "Number of clusters") +
  theme_minimal(); plot_nclust


# plot_nclust_nolegend <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
#   # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
#   geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
#   scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
#   labs(y = "Adjusted rand index (ARI)", x = "Number of clusters") +
#   theme_minimal()+
#   theme(legend.position = "none"); plot_nclust_nolegend


# # save
# png("~/Desktop/plot_nclust.png", width = 14, height = 12, units = 'cm', res = 300)
# plot_nclust
# #don't forget to embed fonts.
# #embed_fonts("plotname.pdf", outfile = "plotname_embed.pdf")
# dev.off()


################################
## Varied Number of Variables ##
################################

## plot netClust versions for different numbers of covariates
res1 <- readRDS("results/results--k_clust-6--n_vars-10--n_bg-5--n_it-20--different--TRUE.rds")
# nice_plot(res1$correct_samples)
nbg_data <- melt(res1$correct_samples[,1:9])
nbg_data$nbg <- 10
res2 <- readRDS("results/results--k_clust-6--n_vars-15--n_bg-5--n_it-20--different--TRUE.rds")
# nice_plot(res2$correct_samples)
nbg_data_temp <- melt(res2$correct_samples[,1:9])
nbg_data_temp$nbg <- 15
nbg_data <- rbind(nbg_data,nbg_data_temp)
res3 <- readRDS("results/results--k_clust-6--n_vars-20--n_bg-5--n_it-20--different--TRUE.rds")
nbg_data_temp <- melt(res3$correct_samples[,1:9])
# nice_plot1 <- nice_plot(res3$correct_samples)
nbg_data_temp$nbg <- 20
nbg_data <- rbind(nbg_data,nbg_data_temp)
res4 <- readRDS("results/results--k_clust-6--n_vars-20--n_bg-5--n_it-20--different--TRUE.rds")
# nice_plot(res4$correct_samples)
nbg_data_temp <- melt(res4$correct_samples[,1:9])
nbg_data_temp$nbg <- 25
nbg_data <- rbind(nbg_data,nbg_data_temp)

colnames(nbg_data)[2] <- "Method"

# remove Mclust (since Gaussian mixture model not interesting for binary data)
nbg_data_reduced <- nbg_data[!nbg_data$Method=="Mclust (cov & var)"&!nbg_data$Method=="Mclust (var)",]
# create labels for publication
levels(nbg_data_reduced$Method) <- c("CANclust (cov. & var.)", "BNMM (cov. & var.)",
                                     "BNMM (var.)", "K-means (cov. & var.)", "K-means (var.)",
                                     "MC", "MC2", "BMM (cov. & var.)", "BMM (var.)")

# choose colors 
# color_list <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")
# cols <- color_list[c(2,9,4,1,10,5,3)]

plot_nvars <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted rand index (ARI)", x = "Number of variables") +
  theme_minimal(); plot_nvars


# plot_nvars_nolegend <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
#   # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
#   geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
#   scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
#   labs(y = "Adjusted rand index (ARI)", x = "Number of variables") +
#   theme_minimal()+
#   theme(legend.position = "none"); plot_nvars_nolegend


##############################
## Varied Number of Samples ##
##############################

## plot netClust versions for different numbers of covariates
res1 <- readRDS("results/results--k_clust-4--n_vars-20--n_bg-5--n_it-600n_samples20.rds")
# nice_plot(res1$correct_samples)
nbg_data <- melt(res1$correct_samples[,1:9])
nbg_data$nbg <- 600
res2 <- readRDS("results/results--k_clust-4--n_vars-20--n_bg-5--n_it-800n_samples20.rds")
# nice_plot(res2$correct_samples)
nbg_data_temp <- melt(res2$correct_samples[,1:9])
nbg_data_temp$nbg <- 800
nbg_data <- rbind(nbg_data,nbg_data_temp)
res3 <- readRDS("results/results--k_clust-4--n_vars-20--n_bg-5--n_it-1000n_samples20.rds")
nbg_data_temp <- melt(res3$correct_samples[,1:9])
# nice_plot1 <- nice_plot(res3$correct_samples)
nbg_data_temp$nbg <- 1000
nbg_data <- rbind(nbg_data,nbg_data_temp)
res4 <- readRDS("results/results--k_clust-4--n_vars-20--n_bg-5--n_it-1200n_samples20.rds")
# nice_plot(res4$correct_samples)
nbg_data_temp <- melt(res4$correct_samples[,1:9])
nbg_data_temp$nbg <- 1200
nbg_data <- rbind(nbg_data,nbg_data_temp)

colnames(nbg_data)[2] <- "Method"

# remove Mclust (since Gaussian mixture model not interesting for binary data)
nbg_data_reduced <- nbg_data[!nbg_data$Method=="Mclust (cov & var)"&!nbg_data$Method=="Mclust (var)",]
# create labels for publication
levels(nbg_data_reduced$Method) <- c("CANclust (cov. & var.)", "BNMM (cov. & var.)",
                                     "BNMM (var.)", "K-means (cov. & var.)", "K-means (var.)",
                                     "MC", "MC2", "BMM (cov. & var.)", "BMM (var.)")

# choose colors 
# color_list <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")
# cols <- color_list[c(2,9,4,1,10,5,3)]

plot_samples <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted rand index (ARI)", x = "Number of samples per cluster") +
  theme_minimal(); plot_samples


# plot_samples_nolegend <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
#   # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
#   geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
#   scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
#   labs(y = "Adjusted rand index (ARI)", x = "Number of samples per cluster") +
#   theme_minimal()+
#   theme(legend.position = "none"); plot_samples_nolegend


####################################
## Combine them ##
####################################

legend_3 <- get_legend(plot_nclust)

plot_test <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted rand index (ARI)", x = "Number of samples per cluster") +
  theme_minimal()+
  theme(legend.position="bottom"); plot_test
legend_4 <- get_legend(plot_test)

# publication figure 4
p_all4 <- ggarrange(plot_nbg, plot_nclust, plot_nvars, plot_samples, nrow = 2, ncol = 2, common.legend = TRUE, legend="bottom"); p_all4

# # save
# png("~/Desktop/plot_all4.png", width = 21, height = 22, units = 'cm', res = 300)
# p_all4
# dev.off()

library("extrafont")
loadfonts()
pdf("~/Desktop/plot_all4.pdf", height = 8, width = 8,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
p_all4
par(op)
dev.off()


# p_all3 <- ggarrange(plot_nbg_nolegend, plot_nclust_nolegend, plot_nclust_nolegend, plot_samples_nolegend, legend_4, nrow = 5); p_all3#, widths=c(0.4,0.4,0.2))
# 
# # save
# png("~/Desktop/plot_all3.png", width = 21, height = 26, units = 'cm', res = 300)
# p_all3
# dev.off()



