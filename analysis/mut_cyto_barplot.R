# plot individual networks
library(ggplot2)
library(ggpubr)
library(extrafont)
library(clustNet)
library(extrafont)

rm(list=ls())
loadfonts()

cluster_results <- readRDS("../results/euler_memberships.rds")
mut_cov_data <- read.table("../data/binary-mutationCovariate-matrix.txt")[,-c(58:61)]

k_clust <- length(cluster_results$DAGs)

# Compute by entropy
entropy <- function(cat.vect){
  # from https://stats.stackexchange.com/questions/221332/variance-of-a-distribution-of-multi-level-categorical-data
  px  = table(cat.vect)/length(cat.vect)
  lpx = log(px, base=2)
  ent = -sum(px*lpx)
  return(ent)
}

# calculate overall entropy
total_entropies <- sapply(mut_cov_data, function(x) entropy(x))
# calculate relative entropy and select variables with the lowest entropy
n_net <- 15 # number of most frequent variables selected variables per cluster
all_var_selected <- c()
cluster_dim <- c()
diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(mut_cov_data))
binary_frequency <- matrix(NA, nrow = k_clust, ncol = ncol(mut_cov_data))
for (nn in 1:k_clust){
  # select data
  cluster_data <- mut_cov_data[cluster_results$clustermembership==nn,]
  cluster_dim[nn] <- dim(cluster_data)[1]
  # calculate frequency
  binary_frequency[nn,] <- sapply(cluster_data, function(x) sum(x))
}

# normalize by cluster size
binary_frequency <- sweep(binary_frequency, 1, cluster_dim, FUN = "/")

mut_cyto_frequ <- binary_frequency[,1:46]
colnames(mut_cyto_frequ) <- colnames(mut_cov_data)[1:46]
# Convert to data frame in long format
binary_frequency_df <- as.data.frame(as.table(mut_cyto_frequ))

colourys<-c("#202020","#774411","#DDAA77","#ed2124","#114477","#CC99BB",
                     "#88CCAA","#117744","#77AADD")
                     
# Update factor levels
levels(binary_frequency_df$Var2)[levels(binary_frequency_df$Var2) == "CG_5"] <- "-5/5q"
levels(binary_frequency_df$Var2)[levels(binary_frequency_df$Var2) == "CG_7"] <- "-7/7q"
levels(binary_frequency_df$Var2)[levels(binary_frequency_df$Var2) == "CG_8"] <- "+8/8q"
levels(binary_frequency_df$Var2)[levels(binary_frequency_df$Var2) == "CG_Y"] <- "-Y"
levels(binary_frequency_df$Var2)[levels(binary_frequency_df$Var2) == "CG_11"] <- "+11/11q"
levels(binary_frequency_df$Var2)[levels(binary_frequency_df$Var2) == "CG_inv16"] <- "inv16"
levels(binary_frequency_df$Var2)[levels(binary_frequency_df$Var2) == "CG_8_21"] <- "t(8;21)"
levels(binary_frequency_df$Var2)[levels(binary_frequency_df$Var2) == "CG_20"] <- "-20/20q"

# Create barplot using ggplot2
p11 <- ggplot(binary_frequency_df, aes(fill=Var1, y=Freq, x=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colourys) +
  xlab("Mutations and cytogenetic alterations") +
  ylab("Frequency") + 
  theme_minimal() +
  labs(fill = "Clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8));p11

# Create barplot using ggplot2
p11 <- ggplot(binary_frequency_df, aes(fill=Var1, y=Freq, x=Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colourys) +
  xlab("Mutations and cytogenetic alterations") +
  ylab("Frequency") + 
  theme_minimal() +
  labs(fill = "Clusters") +
  theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8));p11

cluster_legend <- get_legend(p11)

# save legend
saveRDS(cluster_legend, "../figures/bar_legend.rds")

p_mc <- ggbarplot(binary_frequency_df, "Var2", "Freq",
                fill = "Var1", color = "Var1", palette = colourys,
                label = TRUE, lab.col = NA)+ 
  # label = TRUE, lab.col = "white", lab.pos = "in")+ 
  xlab("Mutations and cytogenetic alterations") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cancer Type"),col = FALSE)+
  theme(legend.position="off",
        axis.line=element_blank(),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        # axis.ticks.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)); p_mc

saveRDS(p_mc, "../figures/mut_cyto_barplot.rds")

pdf("~/Desktop/mut_cyto_barplot.pdf", height = 4, width = 8,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
p11
par(op)
dev.off()


p_mc <- ggbarplot(binary_frequency_df, "Var2", "Freq",
                  fill = "Var1", color = "Var1", palette = colourys,
                  label = TRUE, lab.col = NA)+ 
  # label = TRUE, lab.col = "white", lab.pos = "in")+ 
  xlab("Mutations and cytogenetic alterations") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cancer Type"),col = FALSE)+
  theme(legend.position="right",
        axis.line=element_blank(),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        # axis.ticks.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)); p_mc

saveRDS(p_mc, "../figures/shiny_mut_cyto_barplot.rds")

