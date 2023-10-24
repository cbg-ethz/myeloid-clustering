# Load required libraries
library(networkD3)
library(ggalluvial)
library(extrafont)

# Clear the workspace
rm(list=ls())
loadfonts()

# Read the data
cluster_results <- readRDS("../results/euler_memberships.rds")
sample_memberships <- cluster_results$clustermembership
mutation_covariate_data <- readRDS("../data/risk_scores.rds")

cancer_type <- mutation_covariate_data$Dx

k_clust <- length(cluster_results$DAGs)
sample_memberships <- as.factor(sample_memberships)
levels(sample_memberships) <- LETTERS[1:k_clust]

# Create a data frame with the two classification vectors
data <- data.frame(cancer_type = cancer_type, sample_memberships = sample_memberships)

# create alluvial plot with risk scores
risk_score <- mutation_covariate_data$risk_scores

risk_score <- as.factor(risk_score)
risk_score <- factor(risk_score, levels = c("Adverse (ELN2022)", "Intermediate (ELN2022)", "Favorable (ELN2022)", 
                              "Very High (IPSSM)", "High (IPSSM)", "Moderate High (IPSSM)", 
                              "Moderate Low (IPSSM)", "Low (IPSSM)", "Very Low (IPSSM)"))

# Create a data frame with the two classification vectors
data <- data.frame(cancer_type = cancer_type, sample_memberships = sample_memberships, risk_score = risk_score)

data$cancer_type <- factor(data$cancer_type, levels = c("AML", "MDS", "CMML", "MPN"))

# Create a data frame with the two classification vectors
data <- data.frame(cancer_type = cancer_type, sample_memberships = sample_memberships, risk_score = risk_score)

data <- rbind(data[data$cancer_type=="AML",], data[data$cancer_type=="MDS",])

my_color2 <- c("#DD7788", "#771122")

# # create alluvial plot
# ggplot(data = data,
#        aes(axis1 = cancer_type, axis3 = sample_memberships, axis2 = risk_score)) +
#   geom_alluvium(aes(fill = cancer_type)) +
#   geom_stratum(aes(fill = cancer_type)) +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Survey", "Response"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_manual(values = my_color2) +
#   theme_void() +
#   theme(legend.position = "none") 
# 
# # create alluvial plot
# p1 <- ggplot(data = data,
#        aes(axis1 = cancer_type, axis2 = sample_memberships, axis3 = risk_score)) +
#   geom_alluvium(aes(fill = cancer_type)) +
#   geom_stratum(aes(fill = cancer_type)) +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Survey", "Response"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_manual(values = my_color2) +
#   theme_void() +
#   theme(legend.position = "none"); p1
# 
# 
# library("extrafont")
# loadfonts()
# pdf("~/Desktop/sankey_risk_scores.pdf", height = 7, width = 12,
#     family = "Arial", paper = "special", onefile = FALSE)
# # family = "Times New Roman", paper = "special", onefile = FALSE)
# op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
# p1
# par(op)
# dev.off()

# plot without cluster F
data2 <- data[-which(data$sample_memberships=="F"),]

# create alluvial plot
p2 <- ggplot(data = data2,
             aes(axis1 = cancer_type, axis2 = sample_memberships, axis3 = risk_score)) +
  geom_alluvium(aes(fill = cancer_type)) +
  geom_stratum(aes(fill = cancer_type)) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = my_color2) +
  theme_void() +
  theme(legend.position = "none"); p2

pdf("~/Desktop/sankey_risk_scores_without_G_I.pdf", height = 7, width = 12,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
p2
par(op)
dev.off()
