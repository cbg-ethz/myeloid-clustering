
# Prepare session, load packages
rm(list=ls())
library(survival)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

# define nice colors
colourys <- brewer.pal(7, "RdYlBu")
cols <- colourys[c(1,5,7)]

# load data
output <- read.table("../data/tcga-clinical-information.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

mutCovData <- read.table("../data/binary-mutationCovariate-matrix.txt")

# define useful functions
sinlge_cox_analysis <- function(clusterResMut, mutCovData, output){
  ## Summarize and relabel cluster assignments ##
  cluster.membership <- data.frame("id"=rownames(mutCovData),"groupM"=clusterResMut$clustermembership
                                   ,"groupC"=clusterResCond$clustermembership, "groupP"=clusterResPlain$clustermembership)
  # Change group from 1:22 to A:V
  cluster.membership$groupM <- as.factor(cluster.membership$groupM)
  levels(cluster.membership$groupM) <- LETTERS[1:length(unique(cluster.membership$groupM))]
  
  # Change group from 1:22 to A:V
  cluster.membership$groupC <- as.factor(cluster.membership$groupC)
  levels(cluster.membership$groupC) <- LETTERS[1:length(unique(cluster.membership$groupC))]
  
  # Change group from 1:22 to A:V
  cluster.membership$groupP <- as.factor(cluster.membership$groupP)
  levels(cluster.membership$groupP) <- LETTERS[1:length(unique(cluster.membership$groupP))]
  
  # Merge clinical data with clustering
  clinical <- merge(output,cluster.membership, by="id", sort=FALSE)
  colnames(clinical)[colnames(clinical)=="type"] <- "tissue"
  
  # check how many samples are in each cluster/group
  table(clinical$groupM)
  table(clinical$groupC)
  table(clinical$groupP)
  
  # define continuous and discrete variables
  clinical[,c('age','event','time')] <- lapply(clinical[,c('age','event','time')], as.numeric)
  # Remove subtype from tissue
  clinical[,c('stage','groupM','groupC','groupP','tissue','gender')] <- lapply(clinical[,c('stage','groupM','groupC','groupP','tissue','gender')], factor)
  
  # define empty variable to store results
  coxResults <- matrix(NA,1,4)
  # rownames(coxResults) <- c("Mut","Cond","Plain")
  colnames(coxResults) <- c("LR-corrected","pvalue-corrected","LR","pvalue")
  
  ## cox analysis for all three versions ##
  
  ## mutation only version ##
  
  # Change in likelihood ratio test when variables are added
  # Clinical only
  stageTest <- summary(coxph(Surv(time, event) ~ stage, data = clinical, na.action = "na.omit"))$logtest[1]
  ageTest <- summary(coxph(Surv(time, event) ~ stage + age, data = clinical, na.action = "na.omit"))$logtest[1]
  # Clinical + tissue
  tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + gender, data = clinical, na.action = "na.omit"))$logtest[1]
  # Clinical + tissue + group
  groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + gender + groupM, data = clinical, na.action = "na.omit"))$logtest[1]
  
  LR <- round((groupTest-tissueTest)/2, 1)
  pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$groupM)) - 1, lower.tail = FALSE), 10)
  coxResults[1,1:2] <- c(LR,pvalue)
  paste(round(LR, 1), "&", pvalue, sep = " ")
  
  # Without correction for confounding factors
  
  clusterTest <- summary(coxph(Surv(time, event) ~ groupM, data = clinical, na.action = "na.omit"))$logtest[1]
  LR <- clusterTest/2
  pvalue <- pchisq(q = clusterTest, df = length(levels(clinical$groupM)) - 1, lower.tail = FALSE)
  coxResults[1,3:4] <- c(LR,pvalue)
  paste(round(LR, 1), "&", pvalue, sep = " ")
  
  return(coxResults)
}

new_cox_analysis <- function(clusterResMut, clusterResCond, clusterResPlain, mutCovData, output){
  ## Summarize and relabel cluster assignments ##
  cluster.membership <- data.frame("id"=rownames(mutCovData),"groupM"=clusterResMut$clustermembership
                                   ,"groupC"=clusterResCond$clustermembership, "groupP"=clusterResPlain$clustermembership)
  # Change group from 1:22 to A:V
  cluster.membership$groupM <- as.factor(cluster.membership$groupM)
  levels(cluster.membership$groupM) <- LETTERS[1:length(unique(cluster.membership$groupM))]
  
  # Change group from 1:22 to A:V
  cluster.membership$groupC <- as.factor(cluster.membership$groupC)
  levels(cluster.membership$groupC) <- LETTERS[1:length(unique(cluster.membership$groupC))]
  
  # Change group from 1:22 to A:V
  cluster.membership$groupP <- as.factor(cluster.membership$groupP)
  levels(cluster.membership$groupP) <- LETTERS[1:length(unique(cluster.membership$groupP))]
  
  # Merge clinical data with clustering
  clinical <- merge(output,cluster.membership, by="id", sort=FALSE)
  colnames(clinical)[colnames(clinical)=="type"] <- "tissue"
  
  # check how many samples are in each cluster/group
  table(clinical$groupM)
  table(clinical$groupC)
  table(clinical$groupP)
  
  # define continuous and discrete variables
  clinical[,c('age','event','time')] <- lapply(clinical[,c('age','event','time')], as.numeric)
  # Remove subtype from tissue
  clinical[,c('stage','groupM','groupC','groupP','tissue','gender')] <- lapply(clinical[,c('stage','groupM','groupC','groupP','tissue','gender')], factor)
  
  # define empty variable to store results
  coxResults <- matrix(NA,3,4)
  rownames(coxResults) <- c("Mut","Cond","Plain")
  colnames(coxResults) <- c("LR-corrected","pvalue-corrected","LR","pvalue")
  
  ## cox analysis for all three versions ##
  
  ## mutation only version ##
  
  # Change in likelihood ratio test when variables are added
  # Clinical only
  stageTest <- summary(coxph(Surv(time, event) ~ stage, data = clinical, na.action = "na.omit"))$logtest[1]
  ageTest <- summary(coxph(Surv(time, event) ~ stage + age, data = clinical, na.action = "na.omit"))$logtest[1]
  # Clinical + tissue
  tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + gender, data = clinical, na.action = "na.omit"))$logtest[1]
  # Clinical + tissue + group
  groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + gender + groupM, data = clinical, na.action = "na.omit"))$logtest[1]
  
  LR <- round((groupTest-tissueTest)/2, 1)
  pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$groupM)) - 1, lower.tail = FALSE), 10)
  coxResults[1,1:2] <- c(LR,pvalue)
  paste(round(LR, 1), "&", pvalue, sep = " ")
  
  # Without correction for confounding factors
  
  clusterTest <- summary(coxph(Surv(time, event) ~ groupM, data = clinical, na.action = "na.omit"))$logtest[1]
  LR <- clusterTest/2
  pvalue <- pchisq(q = clusterTest, df = length(levels(clinical$groupM)) - 1, lower.tail = FALSE)
  coxResults[1,3:4] <- c(LR,pvalue)
  paste(round(LR, 1), "&", pvalue, sep = " ")
  
  
  ## conditional covariates version ##
  
  # Change in likelihood ratio test when variables are added
  
  # Clinical only
  stageTest <- summary(coxph(Surv(time, event) ~ stage, data = clinical, na.action = "na.omit"))$logtest[1]
  ageTest <- summary(coxph(Surv(time, event) ~ stage + age, data = clinical, na.action = "na.omit"))$logtest[1]
  # Clinical + tissue
  tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + gender, data = clinical, na.action = "na.omit"))$logtest[1]
  # Clinical + tissue + group
  groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + gender + groupC, data = clinical, na.action = "na.omit"))$logtest[1]
  
  LR <- round((groupTest-tissueTest)/2, 1)
  pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$groupC)) - 1, lower.tail = FALSE), 10)
  coxResults[2,1:2] <- c(LR,pvalue)
  paste(round(LR, 1), "&", pvalue, sep = " ")
  
  # Without correction for confounding factors
  
  clusterTest <- summary(coxph(Surv(time, event) ~ groupC, data = clinical, na.action = "na.omit"))$logtest[1]
  LR <- clusterTest/2
  pvalue <- pchisq(q = clusterTest, df = length(levels(clinical$groupC)) - 1, lower.tail = FALSE)
  coxResults[2,3:4] <- c(LR,pvalue)
  paste(round(LR, 1), "&", pvalue, sep = " ")
  
  ## covariates version ##
  
  # Change in likelihood ratio test when variables are added
  
  # Clinical only
  stageTest <- summary(coxph(Surv(time, event) ~ stage, data = clinical, na.action = "na.omit"))$logtest[1]
  ageTest <- summary(coxph(Surv(time, event) ~ stage + age, data = clinical, na.action = "na.omit"))$logtest[1]
  # Clinical + tissue
  tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + gender, data = clinical, na.action = "na.omit"))$logtest[1]
  # Clinical + tissue + group
  groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + gender + groupP, data = clinical, na.action = "na.omit"))$logtest[1]
  
  LR <- round((groupTest-tissueTest)/2, 1)
  pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$groupP)) - 1, lower.tail = FALSE), 10)
  coxResults[3,1:2] <- c(LR,pvalue)
  paste(round(LR, 1), "&", pvalue, sep = " ")
  
  # Without correction for confounding factors
  
  clusterTest <- summary(coxph(Surv(time, event) ~ groupP, data = clinical, na.action = "na.omit"))$logtest[1]
  LR <- clusterTest/2
  pvalue <- pchisq(q = clusterTest, df = length(levels(clinical$groupP)) - 1, lower.tail = FALSE)
  coxResults[3,3:4] <- c(LR,pvalue)
  paste(round(LR, 1), "&", pvalue, sep = " ")
  
  return(coxResults)
}

 
# Table 1 results
clusterResMut <- readRDS("results/clusterResMut_i200_s1-1.rds")
clusterResCond <- readRDS("results/clusterResCond_i200_s1.rds")
clusterResPlain <- readRDS("results/clusterResPlain_i200_s1-1.rds")
new_cox_analysis(clusterResMut, clusterResCond, clusterResPlain, mutCovData, output)

# Now redo the analysis for multiple seeds

# load the clustering assignments
res1 <- sinlge_cox_analysis(readRDS("results/clusterResCond_i200_s1.rds"), mutCovData, output)
res2 <- sinlge_cox_analysis(readRDS("results/clusterResCond_i200_s2.rds"), mutCovData, output)
res3 <- sinlge_cox_analysis(readRDS("results/clusterResCond_i200_s3.rds"), mutCovData, output)
res4 <- sinlge_cox_analysis(readRDS("results/clusterResCond_i200_s4.rds"), mutCovData, output)
res5 <- sinlge_cox_analysis(readRDS("results/clusterResCond_i200_s5.rds"), mutCovData, output)

res_cond <- rbind(res1, res2, res3, res4, res5)

res1 <- sinlge_cox_analysis(readRDS("results/clusterResMut_i200_s1-1.rds"), mutCovData, output)
res2 <- sinlge_cox_analysis(readRDS("results/clusterResMut_i200_s2-2.rds"), mutCovData, output)
res3 <- sinlge_cox_analysis(readRDS("results/clusterResMut_i200_s3-3.rds"), mutCovData, output)
res4 <- sinlge_cox_analysis(readRDS("results/clusterResMut_i200_s4-4.rds"), mutCovData, output)
res5 <- sinlge_cox_analysis(readRDS("results/clusterResMut_i200_s5-5.rds"), mutCovData, output)

res_mut <- rbind(res1, res2, res3, res4, res5)

res1 <- sinlge_cox_analysis(readRDS("results/clusterResPlain_i200_s1-1.rds"), mutCovData, output)
res2 <- sinlge_cox_analysis(readRDS("results/clusterResPlain_i200_s2-2.rds"), mutCovData, output)
res3 <- sinlge_cox_analysis(readRDS("results/clusterResPlain_i200_s3-3.rds"), mutCovData, output)
res4 <- sinlge_cox_analysis(readRDS("results/clusterResPlain_i200_s4-4.rds"), mutCovData, output)
res5 <- sinlge_cox_analysis(readRDS("results/clusterResPlain_i200_s5-5.rds"), mutCovData, output)

res_plain <- rbind(res1, res2, res3, res4, res5)

# sort and label before plotting
res_all <- cbind(res_cond[,1], res_plain[,1], res_mut[,1])
colnames(res_all) <- c("Covariate-Adjusted", "Variables & Covariates", "Variables")
df_all <- melt(res_all)
colnames(df_all)[2] <- "Data"


# plot to visualize the corrected likelihood ratio for different seeds
temp_plot <- ggplot2::ggplot(data=df_all, aes(x=reorder(Data, -value), y=value, fill = Data, colour = Data))+
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  # geom_point(pch = 21, position = position_jitterdodge(0.15), cex = 0.4, alpha = 0.3) +
  geom_jitter(pch = 21, cex = 0.8, alpha = 0.3) +
  # geom_point(position = position_jitterdodge(), alpha = 0.3, pch = 21, cex = 0.9) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = "Corrected Likelihood Ratio", x = "Method") +
  theme_minimal(); temp_plot

# p_other <- ggarrange(temp_plot, legend_1, widths=c(0.75,0.25))

# save
png("~/Desktop/plot_TCGA_LR.png", width = 12, height = 10, units = 'cm', res = 300) # width = 18, height = 21, units = 'cm', res = 300)
temp_plot
#don't forget to embed fonts.
#embed_fonts("plotname.pdf", outfile = "plotname_embed.pdf")
dev.off()







