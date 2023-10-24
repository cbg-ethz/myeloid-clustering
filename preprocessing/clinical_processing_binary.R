# Prepare session, load packages
rm(list=ls())
setwd("/Users/frbayer/Documents/phd_main/projects/mda_analysis")
library("sjmisc")
library("ggplot2")
require("reshape2")
library("plyr")

# covert AML data to binary data

# read data
mutation_covariate_data <- readRDS("data/aml_data.rds")

# check for NA 
which(mutation_covariate_data==NA) # data is complete

# rename columns
names(mutation_covariate_data)[names(mutation_covariate_data) == 'Gender'] <- 'sex'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'Dx'] <- 'cancer_type'

# get the most frequently mutated genes
mutation_data <- mutation_covariate_data[,-c(1,61:82)]
# save(mutation_data, file = 'data/MDA.rData')

mutation_data_only <- mutation_data[,-c(1)]
n_genes <- ncol(mutation_data_only)
mutation_frequency <- c()
for (i in 1:n_genes){
  mutation_frequency[i] <- table(mutation_data_only[,i])[2]
}
mutation_frequency[which(is.na(mutation_frequency))] <- 0
limit <- 13 # minimum of mutations in a gene to be considered
excluded_genes <- setdiff((1:length(mutation_frequency)),which(mutation_frequency>limit))

# number of genes that are considred
sum(mutation_frequency>limit)
# names of considered genes
colnames(mutation_data_only)[which(mutation_frequency>limit)]

mutation_data <- mutation_data[,c(1,which(mutation_frequency>limit)+1)]

save(mutation_data, file = 'data/MDA.rData')

# dichotomize covariates
temp_mutation_covariate_data <- mutation_covariate_data[,-c(1,2, excluded_genes+2)]

# create continuous copy for plotting
continuous_mutation_covariate_data <- mutation_covariate_data[,-c(1,2, excluded_genes+2)]

# age
# temp_mutation_covariate_data$age <- dicho(mutation_covariate_data$age)

# age by cancer type
age_medians <- ddply(mutation_covariate_data, "cancer_type", summarise, grp.mean=median(age))
head(age_medians)
ggplot(mutation_covariate_data, aes(x=age, color=cancer_type)) +
  geom_density() +
  # scale_x_continuous(limits=c(0,100))+
  geom_vline(data=age_medians, aes(xintercept=grp.mean, color=cancer_type),
             linetype="dashed")

# age by cancer type
age_medians <- ddply(mutation_covariate_data, "sex", summarise, grp.mean=median(age))
head(age_medians)
ggplot(mutation_covariate_data, aes(x=age, color=sex)) +
  geom_density() +
  # scale_x_continuous(limits=c(0,100))+
  geom_vline(data=age_medians, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")

# age binarization
temp_mutation_covariate_data$age[mutation_covariate_data$cancer_type=="AML"] <- 
  dicho(mutation_covariate_data$age[mutation_covariate_data$cancer_type=="AML"])
temp_mutation_covariate_data$age[mutation_covariate_data$cancer_type=="MDS"] <- 
  dicho(mutation_covariate_data$age[mutation_covariate_data$cancer_type=="MDS"])
temp_mutation_covariate_data$age[mutation_covariate_data$cancer_type=="MPN"] <- 
  dicho(mutation_covariate_data$age[mutation_covariate_data$cancer_type=="MPN"])
temp_mutation_covariate_data$age[mutation_covariate_data$cancer_type=="CMML"] <- 
  dicho(mutation_covariate_data$age[mutation_covariate_data$cancer_type=="CMML"])
temp_mutation_covariate_data$age <- temp_mutation_covariate_data$age-1

# blood
temp_mutation_covariate_data$PB_ANC <- dicho(mutation_covariate_data$PB_ANC)
temp_mutation_covariate_data$PB_AMC <- dicho(mutation_covariate_data$PB_AMC)
temp_mutation_covariate_data$PB_ALC <- dicho(mutation_covariate_data$PB_ALC)
temp_mutation_covariate_data$PB_WBC <- dicho(mutation_covariate_data$PB_WBC)
temp_mutation_covariate_data$PB_PLT <- dicho(mutation_covariate_data$PB_PLT)
temp_mutation_covariate_data$PB_HGB <- dicho(mutation_covariate_data$PB_HGB)
temp_mutation_covariate_data$PB_Blast <- as.numeric(cut(mutation_covariate_data$BM_Blast, breaks = c(-1,9,19,100)))-1

# bone
temp_mutation_covariate_data$BM_Blast <- dicho(mutation_covariate_data$BM_Blast)
temp_mutation_covariate_data$BM_Mono <- dicho(mutation_covariate_data$BM_Mono)
temp_mutation_covariate_data$BM_Gran <- dicho(mutation_covariate_data$BM_Gran)

# create binary columns for cancer type
cancer_type <- matrix(0,dim(temp_mutation_covariate_data)[1],4)
colnames(cancer_type) <- c("AML", "CMML", "MDS", "MPN")
cancer_type[,1][which(mutation_covariate_data$cancer_type=="AML")] <- 1
cancer_type[,2][which(mutation_covariate_data$cancer_type=="CMML")] <- 1
cancer_type[,3][which(mutation_covariate_data$cancer_type=="MDS")] <- 1
cancer_type[,4][which(mutation_covariate_data$cancer_type=="MPN")] <- 1

temp_mutation_covariate_data <- cbind(temp_mutation_covariate_data, cancer_type)

# put age in the last col of the matrix
index_age <- which(colnames(temp_mutation_covariate_data)=="age")
temp_mutation_covariate_data <- temp_mutation_covariate_data[,c(1:(index_age-1),(index_age+1):ncol(temp_mutation_covariate_data), index_age)]
# temp_mutation_covariate_data[ , c(index_age,ncol(temp_mutation_covariate_data))]  <- temp_mutation_covariate_data[ , c(ncol(temp_mutation_covariate_data),index_age)]
# colnames(temp_mutation_covariate_data)[c(index_age,ncol(temp_mutation_covariate_data))]  <- colnames(temp_mutation_covariate_data)[c(ncol(temp_mutation_covariate_data),index_age)]

# sex (male=1)
temp_mutation_covariate_data$sex <- factor(mutation_covariate_data$sex, levels=c("Female", "Male"), labels=c(0, 1))

# # rearrange
# temp_mutation_covariate_data[,c(1,3:60,62:81,2,82,61)]
# temp_mutation_covariate_data[,c(1,3:60,62:81,2,82,61)]

# remove PB_AlC, as not really informative
temp_mutation_covariate_data$PB_ALC <- NULL

# put causal covariates at the end of the table
temp_mutation_covariate_data <- temp_mutation_covariate_data[,c(1:57,59:63,58)]

# change order of PB and CG values
temp_mutation_covariate_data <- temp_mutation_covariate_data[,c(1:38,45:52,39:44,53:63)]

# save results
write.table(temp_mutation_covariate_data, file="data/binary-mutationCovariate-matrix.txt", quote=FALSE, sep="\t")

# investigate the binarization 

# get blood colnames
colnames(continuous_mutation_covariate_data)[c(60:66)]
blood_values <- continuous_mutation_covariate_data[,c(60:66)]
blood_values_df <- melt(blood_values)

blood_medians <- ddply(blood_values_df, "variable", summarise, grp.mean=median(value))
head(blood_medians)

# add median lines
ggplot(blood_values_df, aes(x=value, color=variable)) +
  geom_density() +
  scale_x_continuous(limits=c(0,10))+
  geom_vline(data=blood_medians, aes(xintercept=grp.mean, color=variable),
             linetype="dashed")

# get bone colnames
colnames(continuous_mutation_covariate_data)[75:79]
bone_values <- continuous_mutation_covariate_data[,c(75:79)]
bone_values_df <- melt(bone_values)

bone_medians <- ddply(bone_values_df, "variable", summarise, grp.mean=median(value))
head(bone_medians)

# add median lines
ggplot(bone_values_df, aes(x=value, color=variable)) +
  geom_density() +
  scale_x_continuous(limits=c(0,10))+
  geom_vline(data=bone_medians, aes(xintercept=grp.mean, color=variable),
             linetype="dashed")

# get blood and bone colnames
colnames(continuous_mutation_covariate_data)[c(66,75)]
bb_values <- continuous_mutation_covariate_data[,c(66,75)]
bb_values_df <- melt(bb_values)

bb_medians <- ddply(bb_values_df, "variable", summarise, grp.mean=median(value))
head(bb_medians)

# add median lines
ggplot(bb_values_df, aes(x=value, color=variable)) +
  geom_density() +
  scale_x_continuous(limits=c(0,100))+
  geom_vline(data=bb_medians, aes(xintercept=grp.mean, color=variable),
             linetype="dashed") +
  geom_vline(aes(xintercept=10),
             linetype="dashed") +
  geom_vline(aes(xintercept=20),
           linetype="dashed")

