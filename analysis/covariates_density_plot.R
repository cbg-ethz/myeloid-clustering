# Load required libraries
library("sjmisc")
library("ggplot2")
require("reshape2")
library("plyr")
library("ggpubr")
library("extrafont")

# Clear the workspace
rm(list=ls())
loadfonts()

# read data
mutation_covariate_data <- readRDS("../data/aml_data.rds")

# check for NA 
# which(mutation_covariate_data==NA) # data is complete

names(mutation_covariate_data)[names(mutation_covariate_data) == 'Gender'] <- 'sex'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'Dx'] <- 'cancer_type'

# dichotomize covariates
temp_mutation_covariate_data <- mutation_covariate_data[,-c(1,2)]

# create continuous copy for plotting
continuous_mutation_covariate_data <- mutation_covariate_data[,-c(1,2)]

# age
# temp_mutation_covariate_data$age <- dicho(mutation_covariate_data$age)

# age by cancer type
age_medians <- ddply(mutation_covariate_data, "cancer_type", summarise, grp.mean=median(age))
# head(age_medians)
age_medians_table <- age_medians
colnames(age_medians_table) <- c("Cancer type","Median")
knitr::kable(t(age_medians_table), caption = "Median age per cancer type")
p_age <- ggplot(mutation_covariate_data, aes(x=age, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
             linetype="dashed") +
  theme_minimal() +
  labs(color='Cancer Type') +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  ggtitle("Age Distribution")+
  ylab("Density")+
  xlab("Age");p_age

blood_values <- mutation_covariate_data[,c(2,68)]
blood_values_df <- melt(blood_values)

p_pb <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  scale_x_continuous(limits=c(0,10))+
  theme_minimal() +
  ggtitle(as.character(blood_values_df[1,2])) +
  labs(color='Cancer Type') +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+ 
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  ggtitle("PB Blasts")+
  ylab("Density")+
  xlab("Value");p_pb

blood_values <- mutation_covariate_data[,c(2,67)]
blood_values_df <- melt(blood_values)

p2 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
             linetype="dashed") +
  geom_polygon(data = data.frame(x = c(11.9, 16.9, 16.9, 11.9), y = c(-Inf, -Inf, Inf, Inf)), aes(x = x, y = y), color = NA, fill = "red", alpha = 0.1) +
  scale_x_continuous(limits=c(0,20))+
  theme_minimal() +
  labs(color='Cancer Type') +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+ 
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  ggtitle("Hemoglobin")+
  ylab("Density")+
  xlab("Value"); p2

blood_values <- mutation_covariate_data[,c(2,66)]
blood_values_df <- melt(blood_values)

p3 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  geom_polygon(data = data.frame(x = c(152, 361, 361, 152), y = c(-Inf, -Inf, Inf, Inf)), aes(x = x, y = y), color = NA, fill = "red", alpha = 0.1) +
  scale_x_continuous(limits=c(0,500))+
  theme_minimal() +
  labs(color='Cancer Type') +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  ggtitle("Platelets")+
  ylab("Density")+
  xlab("Value"); p3

blood_values <- mutation_covariate_data[,c(2,65)]
blood_values_df <- melt(blood_values)

# blood_values <- mutation_covariate_data[,c(2,68)]
# blood_values_df <- melt(blood_values)

p4 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  geom_polygon(data = data.frame(x = c(4.5, 11.5, 11.5, 4.5), y = c(-Inf, -Inf, Inf, Inf)), aes(x = x, y = y), color = NA, fill = "red", alpha = 0.1) +
  scale_x_continuous(limits=c(0,30))+
  theme_minimal() +
  labs(color='Cancer Type') +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  ggtitle("White Blood Cells")+
  ylab("Density")+
  xlab("Value");p4

blood_values <- mutation_covariate_data[,c(2,64)]
blood_values_df <- melt(blood_values)

p5 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  scale_x_continuous(limits=c(0,10))+
  theme_minimal() +
  labs(color='Cancer Type') +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  geom_polygon(data = data.frame(x = c(1, 4.8, 4.8, 1), y = c(-Inf, -Inf, Inf, Inf)), aes(x = x, y = y), color = NA, fill = "red", alpha = 0.1) +
  ggtitle("Absolute Lymphocyte Count")+
  ylab("Density")+
  xlab("Value"); p5


blood_values <- mutation_covariate_data[,c(2,63)]
blood_values_df <- melt(blood_values)

p6 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  scale_x_continuous(limits=c(0,10))+
  theme_minimal() +
  labs(color='Cancer Type') +
  ggtitle(as.character(blood_values_df[1,2])) +
  geom_polygon(data = data.frame(x = c(2, 8, 8, 2), y = c(-Inf, -Inf, Inf, Inf)), aes(x = x, y = y), color = NA, fill = "red", alpha = 0.1) +
  # scale_colour_brewer(type = "seq", palette = "Spectral") +
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  ggtitle("Absolute Monocyte Count")+
  ylab("Density")+
  xlab("Value"); p6


blood_values <- mutation_covariate_data[,c(2,62)]
blood_values_df <- melt(blood_values)

p7 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  geom_polygon(data = data.frame(x = c(2.5, 6, 6, 2.5), y = c(-Inf, -Inf, Inf, Inf)), aes(x = x, y = y), color = NA, fill = "red", alpha = 0.1) +
  scale_x_continuous(limits=c(0,20)) +
  theme_minimal() +
  labs(color='Cancer Type') +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  ggtitle("Absolute Neutrophil Count")+
  ylab("Density")+
  xlab("Value");p7 


# ggarrange(p2,p3,p4,p5,p6,p7, nrow = 2)
ggarrange(p3,p2,p4, p7, p5, p6, common.legend = TRUE, legend="bottom")

pdf("~/Desktop/blood_values.pdf", height = 6, width = 10,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
ggarrange(p3,p2,p4, p7, p5, p6, common.legend = TRUE, legend="bottom")
par(op)
dev.off()

p_blood <- ggarrange(p3, p2, p4, p7, p5, p6, common.legend = TRUE, legend="bottom")


## --------------------------------------------------


blood_values <- mutation_covariate_data[,c(2,77)]
blood_values_df <- melt(blood_values)

p_bm <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(adjust = 1.6) +
  # scale_x_continuous(limits=c(0,100))+
  geom_vline(data=age_medians, aes(xintercept=10), linetype="dashed") +
  geom_vline(data=age_medians, aes(xintercept=20), linetype="dashed") +
  scale_x_continuous(limits=c(0,100))+
  theme_minimal() +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  ggtitle("BM Blast")+
  labs(color='Cancer Type') +
  # geom_vline(data=bb_third1[2,], aes(xintercept=grp.mean),colour="#BB0000",
  #            linetype="dashed") +
  # geom_vline(data=bb_third2[2,], aes(xintercept=grp.mean),colour="#BB0000",
  #            linetype="dashed") +
  ylab("Density")+
  xlab("Value"); p_bm

p_bm2 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(adjust = 1.6) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=10), linetype="dashed") +
  # geom_vline(data=age_medians, aes(xintercept=20), linetype="dashed") +
  scale_x_continuous(limits=c(0,100))+
  theme_minimal() +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  ggtitle("BM Blast")+
  labs(color='Cancer Type') +
  # geom_vline(data=bb_third1[2,], aes(xintercept=grp.mean),colour="#BB0000",
  # linetype="dashed") +
  # geom_vline(data=bb_third2[2,], aes(xintercept=grp.mean),colour="#BB0000",
  # linetype="dashed") +
  ylab("Density")+
  xlab("Value"); p_bm2

ggarrange(p_pb,p_bm2, common.legend = TRUE, legend="bottom")
ggarrange(p_pb,p_bm, common.legend = TRUE, legend="bottom")

# p_blasts <- ggarrange(p_pb,p_bm2, common.legend = TRUE, legend="bottom", ncol = 1)
p_blasts <- ggarrange(p_pb,p_bm2, common.legend = TRUE, legend="none", ncol = 1)

p_blood2 <- ggarrange(p_pb, p3, p2, p4, p7, p6, common.legend = TRUE, legend="bottom")

# ---------------------------------------------

blood_values <- mutation_covariate_data[,c(2,78)]
blood_values_df <- melt(blood_values)

p2 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  scale_x_continuous(limits=c(0,30))+
  # geom_polygon(data = data.frame(x = c(2, 8, 8, 2), y = c(-Inf, -Inf, Inf, Inf)), aes(x = x, y = y), color = NA, fill = "red", alpha = 0.1) +
  theme_minimal() +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  labs(color='Cancer Type') +
  ggtitle("Monocytes in the Bone Marrow")+
  ylab("Density")+
  xlab("Value"); p2

blood_values <- mutation_covariate_data[,c(2,79)]
blood_values_df <- melt(blood_values)

p3 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  scale_x_continuous(limits=c(0,100))+
  theme_minimal() +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  labs(color='Cancer Type') +
  ggtitle("Bone Marrow Granulocytes")+
  ylab("Density")+
  xlab("Value"); p3


blood_values <- mutation_covariate_data[,c(2,80)]
blood_values_df <- melt(blood_values)

p4 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  scale_x_continuous(limits=c(0,1))+
  theme_minimal() +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  labs(color='Cancer Type') +
  ggtitle("Bone Marrow Auer Rods")+
  ylab("Density")+
  xlab("Value"); p4

blood_values <- mutation_covariate_data[,c(2,81)]
blood_values_df <- melt(blood_values)

p5 <- ggplot(blood_values_df, aes(x=value, color=cancer_type)) +
  geom_density(key_glyph = draw_key_path) +
  # scale_x_continuous(limits=c(0,100))+
  # geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
  # linetype="dashed") +
  scale_x_continuous(limits=c(0,1))+
  theme_minimal() +
  ggtitle(as.character(blood_values_df[1,2])) +
  # scale_colour_brewer(type = "seq", palette = "Spectral")+
  scale_color_manual(values=c("#DD7788", "#117777", "#771122", "#DDDD77"))+
  labs(color='Cancer Type') +
  ggtitle("Bone Marrow Dysplasia")+
  ylab("Density")+
  xlab("Value"); p5

# p1

ggarrange(p2,p3,p4,p5, nrow = 2, ncol = 2)
p_bone <- ggarrange(p2,p3,p4,p5, nrow = 2, ncol = 2, common.legend = TRUE, legend="none")

p_bone2 <- ggarrange(p_bm2, p2,p3, nrow = 1, ncol = 3, common.legend = TRUE, legend="none")

# p_all_vars <- ggarrange(ggarrange(p_blasts, p_bone, ncol = 2, labels = c("a", "b"),widths = c(1,2)),p_blood, # Second row with box and dot plots
#                         nrow = 2, 
#                         labels = c(NA,"c"),                                        # Labels of the scatter plot
#                         common.legend = TRUE, legend="bottom");p_all_vars
# 
# 
# p_all_vars <- ggarrange(ggarrange(p_blasts+theme(plot.background = element_rect(color = "black")), p_bone+ theme(plot.background = element_rect(color = "black")), ncol = 2, labels = c("a", "b"),widths = c(1,2)),p_blood+theme(plot.background = element_rect(color = "black")), # Second row with box and dot plots
#                         nrow = 2, 
#                         labels = c(NA,"c"),                                        # Labels of the scatter plot
#                         common.legend = TRUE, legend="bottom");p_all_vars
# 
# pdf("~/Desktop/cov_values.pdf", height = 8, width = 10,
#     family = "Arial", paper = "special", onefile = FALSE)
# # family = "Times New Roman", paper = "special", onefile = FALSE)
# op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
# p_all_vars
# par(op)
# dev.off()

p_all_vars2 <- ggarrange(p_bone2,p_blood2, # Second row with box and dot plots
                        nrow = 2, heights = c(1,2),
                        labels = c("a","b"),                                        # Labels of the scatter plot
                        common.legend = TRUE, legend="bottom");p_all_vars2

pdf("~/Desktop/cov_values2.pdf", height = 7, width = 10,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
p_all_vars2
par(op)
dev.off()


