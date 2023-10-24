# Print summary tables of the survival analysis

# Load required libraries
library(survival)
library(survminer)
library(dplyr)
library(xtable)

# Clear the workspace
rm(list=ls())

# read data
clinical <- readRDS("../data/os_data.rds")
cluster_results <- readRDS("../results/euler_memberships_binary_edgepf4.rds")
mutation_covariate_data <- readRDS("../data/aml_data.rds")
risk_scores <- readRDS("../data/risk_scores.rds")

# merge features
clinical$group <- as.factor(cluster_results$clustermembership)
levels(clinical$group) <- LETTERS[1:9]

clinical$type <- mutation_covariate_data$Dx
clinical$gender <- mutation_covariate_data$Gender
clinical$age <- mutation_covariate_data$age
clinical$risk_scores <- risk_scores$risk_scores

# print the summary tobles for the survival analysis by cluster
surv_table <- summary(survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical))$table[,-c(2,3)]

clust_name <- unname(sapply(rownames(surv_table), function(x) substr(x,7,7))) 
surv_table <- data.frame(surv_table)
surv_table <- cbind(clust_name,surv_table)

print(xtable(surv_table, type = "latex", caption="Summary table of survival analysis by each cluster.", digits=1),include.rownames=FALSE)

# print the summary tobles for the survival analysis by cluster and cancer type
surv_table <- summary(survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group+type, data = clinical))$table[,-c(2,3)]

clust_name <- unname(sapply(rownames(surv_table), function(x) substr(x,7,7))) 
ct_name <- unname(sapply(rownames(surv_table), function(x) substr(x,15,18))) 
surv_table <- data.frame(surv_table)
surv_table <- cbind(clust_name,ct_name,surv_table)

print(xtable(surv_table, type = "latex", caption="Summary table of survival analysis by each cluster and cancer type.", digits=1),include.rownames=FALSE)



# print the summary tobles for the survival analysis by risk score and cancer type
surv_table <- summary(survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ risk_scores+type, data = clinical))$table[,-c(2,3)]

clust_name <- unname(sapply(rownames(surv_table), function(x) substr(x,13,34)))
ct_name <- unname(sapply(rownames(surv_table), function(x) substr(x,31,45)))
surv_table <- data.frame(surv_table)
surv_table <- cbind(clust_name,ct_name,surv_table)

print(xtable(surv_table, type = "latex", caption="Summary table of survival analysis by each cluster and cancer type.", digits=1),include.rownames=FALSE)


