
# Prepare session, load packages
rm(list=ls())
library(survival)

## load data ##
output <- read.table("data/tcga-clinical-information.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

BNMM_var_memberships <- readRDS("results/BNMM_var_memberships.rds")
adjust_memberships <- readRDS("results/adjust_memberships.rds")
BNMM_var_co_memberships <- readRDS("results/BNMM_var-co_memberships.rds")

mutCovData <- read.table("data/binary-mutationCovariate-matrix.txt")

## Summarize and relabel cluster assignments ##
cluster.membership <- data.frame("id"=rownames(mutCovData),"groupM"=BNMM_var_memberships$clustermembership
                                 ,"groupC"=adjust_memberships$clustermembership, "groupP"=BNMM_var_co_memberships$clustermembership)
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

# display summary
coxResults

