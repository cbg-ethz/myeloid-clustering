
# Prepare session, load packages
rm(list=ls())
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(clustNet)
library(reshape2)
# library(extrafont)

###########################
## All Methods displayed ##
###########################

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
color_list <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")
cols <- color_list[c(2,9,4,10,1,7,5,8,3)]

# publication figure 2

# remove Mclust (since Gaussian mixture model not interesting for binary data)
nbg_data_reduced <- nbg_data[!nbg_data$Method=="Mclust (cov & var)"&!nbg_data$Method=="Mclust (var)",]
# create labels for publication
levels(nbg_data_reduced$Method) <- c("CANclust (mut. & cov.)", "BNMM (mut. & cov.)",
                             "BNMM (mut.)", "K-means (mut. & cov.)", "K-means (mut.)",
                             "MC", "MC2", "BMM (mut. & cov.)", "BMM (mut.)")

color_list <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")
color_list2 <- RColorBrewer::brewer.pal(n = 7, name = "BrBG")
cols <- color_list[c(2,9,10,3,4,6,7)]
cols[c(4,5)] <- color_list2[c(5,6)]

# color_list <- RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")
color_list <- RColorBrewer::brewer.pal(n = 10, name = "Paired")
cols <- color_list[c(6,1,2,7,8,3,4)]
color_list2 <- RColorBrewer::brewer.pal(n = 6, name = "Greens")
color_list2 <- RColorBrewer::brewer.pal(n = 9, name = "Greens")
cols[c(6,7)] <- color_list2[c(3,5)]
cols[c(6,7)] <- color_list2[c(7,9)]

# cols <- c("#E31A1C", "#114477", "#77AADD", "#771155", "#CC99BB", "#117744", "#88CCAA")
# cols <- c("#E31A1C", "#114477", "#77AADD", "#771122","#DD7788", "#117777","#77CCCC")
# cols <- RColorBrewer::brewer.pal(n = 8, name = "RdBu")
# colourys<-c("#202020","#771122","#AA4455","#DD7788","#774411","#AA7744",
#                      "#DDAA77","#777711","#AAAA44","#DDDD77","#117744","#44AA77",
#                      "#88CCAA","#117777","#44AAAA","#77CCCC","#114477","#4477AA",
#                      "#77AADD","#771155","#AA4488","#CC99BB")
cols <- c("#E31A1C", "#117777", "#77CCCC", "#771122","#DD7788", "#777711", "#DDDD77")
# colourys2 <- c("#DD7788", "#771122", "#117777", "#DDDD77")
# colourys<-c("#202020","#771122","#AA4455","#DD7788","#774411","#AA7744",
#                      "#DDAA77","#777711","#AAAA44","#DDDD77","#117744","#44AA77",
                     # "#88CCAA","#117777","#44AAAA","#77CCCC","#114477","#4477AA",
                     # "#77AADD","#771155","#AA4488","#CC99BB")
                     

plot_nbg <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Number of Covariates") +
  theme_minimal(); plot_nbg

plot_nbg_nolegend <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Adjusted Rand Index (ARI)", x = "Number of Covariates") +
  theme_minimal()+
  theme(legend.position = "none"); plot_nbg_nolegend


# # save
# png("~/Desktop/presentation_plot_nbg.png", width = 18, height = 7, units = 'cm', res = 300)
# plot_nbg
# dev.off()

library("extrafont")
loadfonts()
pdf("~/Desktop/plot_nbg.pdf", height = 3.5, width = 7.5,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
plot_nbg
par(op)
dev.off()


# library("ggsci")
# plot_nbg <- ggplot2::ggplot(data=nbg_data_reduced, aes(x=factor(nbg), y=value, fill = Method, colour = Method))+
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
#   # geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
#   geom_point(position = position_jitterdodge(), alpha = 0.5, pch = 21, cex = 0.9) +
#   scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
#   labs(y = "Adjusted Rand Index (ARI)", x = "Number of Covariates") +
#   theme_minimal()+scale_color_lancet()+scale_fill_lancet(); plot_nbg



