# get bar plots of cluster assignment

# Load required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Clear the workspace
rm(list=ls())

# read data
cluster_results <- readRDS("../results/euler_memberships.rds")
mutation_covariate_data <- readRDS("../data/aml_data.rds")

table(cluster_results$clustermembership)

k_clust <- length(levels(as.factor(cluster_results$clustermembership)))

cancer_table <- matrix(NA, k_clust, 4)

colnames(cancer_table) <- c("AML", "MDS", "MPN", "CMML")
rownames(cancer_table) <- LETTERS[1:k_clust]

for (ii in 1:k_clust){
  count_aml <- length(which(mutation_covariate_data$Dx[which(cluster_results$clustermembership==ii)]=="AML"))
  count_mds <- length(which(mutation_covariate_data$Dx[which(cluster_results$clustermembership==ii)]=="MDS"))
  count_mpn <- length(which(mutation_covariate_data$Dx[which(cluster_results$clustermembership==ii)]=="MPN"))
  count_cmml <- length(which(mutation_covariate_data$Dx[which(cluster_results$clustermembership==ii)]=="CMML"))
  
  cancer_table[ii,] <- c(count_aml, count_mds, count_mpn, count_cmml)
}


mycolor <- c("#DD7788", "#771122", "#DDDD77", "#117777")

p2 <- ggbarplot(melt(cancer_table), "Var1", "value",
          fill = "Var2", color = "Var2", palette = mycolor,
          label = TRUE, lab.col = NA)+ 
  # label = TRUE, lab.col = "white", lab.pos = "in")+ 
  xlab("Clusters") + 
  # ylab("Number of Patients") +
  ylab("Patients per cluster") +
  guides(fill=guide_legend(title="Cancer Type"),col = FALSE)+
  theme(legend.position="off",
        axis.line=element_blank(),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank()); p2

# save plot
saveRDS(p2, "../figures/barplot.rds")

p22 <- ggbarplot(melt(cancer_table), "Var1", "value",
                 fill = "Var2", color = "Var2", palette = mycolor,
                 label = TRUE, lab.col = NA)+ 
  # label = TRUE, lab.col = "white", lab.pos = "in")+ 
  xlab("Clusters") + 
  # ylab("Number of Patients") +
  ylab("") +
  guides(fill=guide_legend(title="Cancer type"),col = FALSE)+
  theme(legend.position="bottom",
        axis.line=element_blank(),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank()); p22

cluster_legend <- get_legend(p22)

# save legend
saveRDS(cluster_legend, "../figures/bar_legend_ct.rds")


# df_cancer_table <- melt(cancer_table)
# 
# df_cancer_table$Var1 <- factor(df_cancer_table$Var1,                                    # Change ordering manually
#                                levels = c("A","H","F","I","C","E","D","B","G"))
# 
# p2 <- ggbarplot(df_cancer_table, "Var1", "value",
#                 fill = "Var2", color = "Var2", palette = mycolor,
#                 label = TRUE, lab.col = NA)+ 
#   # label = TRUE, lab.col = "white", lab.pos = "in")+ 
#   xlab("Clusters") + 
#   # ylab("Number of Patients") +
#   ylab("") +
#   guides(fill=guide_legend(title="Cancer Type"),col = FALSE)+
#   theme(legend.position="off",
#         axis.line=element_blank(),
#         axis.text.y=element_blank(),  #remove y axis labels
#         axis.ticks.y=element_blank(),
#         axis.ticks.x=element_blank()); p2
# 
# # save plot
# saveRDS(p2, "../figures/barplot.rds")


