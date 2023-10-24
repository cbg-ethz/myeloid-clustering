# Survival analysis

# Load required libraries
library(survival)
library(RColorBrewer)
library(survminer)
library(dplyr)

# Clear the workspace
rm(list=ls())

# read data
clinical <- readRDS("../data/os_data.rds")
cluster_results <- readRDS("../results/euler_memberships.rds")
mutation_covariate_data <- readRDS("../data/aml_data.rds")

# merge features
clinical$group <- as.factor(cluster_results$clustermembership)
levels(clinical$group) <- LETTERS[1:9]

clinical$type <- mutation_covariate_data$Dx
clinical$gender <- mutation_covariate_data$Gender
clinical$age <- mutation_covariate_data$age


# Change in likelihood ratio test when variables are added
# Clinical + tissue
tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ type + age, data = clinical, na.action = "na.omit"))$logtest[1]
# Clinical + tissue + group
groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ group + type + age, data = clinical, na.action = "na.omit"))$logtest[1]

LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$group)) - 1, lower.tail = FALSE), 10)
# coxResults[1,1:2] <- c(LR,pvalue)
paste(round(LR, 1), "&", pvalue, sep = " ")


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

# Kaplan-Meier curve for groups
colourysdots <- c("#202020","#774411","#DDAA77","#ed2124","#114477","#CC99BB",
                               "#88CCAA","#117744","#77AADD")

p_surv_temp <- ggsurvplot(survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical), clinical, 
                          palette = colourysdots,
                          legend = "none", 
                          xlab="Time (years)", 
                          xlim = c(0,10),         # present narrower X axis, but not affect
                          size=0.6)

p_surv <- p_surv_temp$plot +scale_y_continuous(
  breaks = c(0,0.2,0.4,0.6,0.8,1),
  expand = c(0.01, 0)
)+scale_x_continuous(
  breaks = c(0,2,4,6,8,10),
  expand = c(0.01, 0)
);p_surv

# save plot
saveRDS(p_surv, "../figures/km_plot.rds")

os <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical)
km_plot <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = colourysdots, # personalized colours
  legend.labs = sort(as.character(unique(clinical$group))), 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
); km_plot

cluster_legend <- get_legend(km_plot$plot)

# save legend
saveRDS(cluster_legend, "../figures/km_legend.rds")

os <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ type, data = clinical)
km_plot <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = c("#DD7788", "#117777", "#771122", "#DDDD77"), # personalized colours
  legend.labs = sort(as.character(unique(clinical$type))), 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
); km_plot

cluster_legend <- get_legend(km_plot$plot)

# save legend
saveRDS(cluster_legend, "../figures/km_legend_ct.rds")


# Change in likelihood ratio test when variables are added
# Clinical only
stageTest <- summary(coxph(Surv(time, event) ~ group, data = clinical, na.action = "na.omit"))$logtest[1]
ageTest <- summary(coxph(Surv(time, event) ~ group + age, data = clinical, na.action = "na.omit"))$logtest[1]
# Clinical + tissue
tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ type + age + gender, data = clinical, na.action = "na.omit"))$logtest[1]
# Clinical + tissue + group
groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ group + type + age + gender, data = clinical, na.action = "na.omit"))$logtest[1]

LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$group)) - 1, lower.tail = FALSE), 10)
# coxResults[1,1:2] <- c(LR,pvalue)
paste(round(LR, 1), "&", pvalue, sep = " ")

# investigate subset of patients diagnosed with AML and MDS 
risk_scores <- readRDS("../data/risk_scores.rds")

clinical$risk_scores <- risk_scores$risk_scores

clinical_aml_mds <- as.data.frame(clinical)

clinical_aml_mds <- clinical_aml_mds[!is.na(clinical$risk_scores),]

# check which classification is more predictive in survival
# check IPSSR_ELN vs our clustering, given the WHO_2016
summary(coxph(Surv(time, as.numeric(event)) ~ group + type + age + gender, data = clinical_aml_mds, na.action = "na.omit"))$logtest[1]
summary(coxph(Surv(time, as.numeric(event)) ~ risk_scores + type + age + gender, data = clinical_aml_mds, na.action = "na.omit"))$logtest[1]


# Clinical + tissue
tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ risk_scores + type + age + gender, data = clinical_aml_mds, na.action = "na.omit"))$logtest[1]
# Clinical + tissue + group
groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ group + risk_scores + type + age + gender, data = clinical_aml_mds, na.action = "na.omit"))$logtest[1]

LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$group)) - 1, lower.tail = FALSE), 10)
# coxResults[1,1:2] <- c(LR,pvalue)
paste(round(LR, 1), "&", pvalue, sep = " ")

