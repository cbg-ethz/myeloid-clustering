# Survival analysis by cancer type

# Load required libraries
library(survival)
library(RColorBrewer)
library(survminer)
library(dplyr)

# Clear the workspace
rm(list=ls())

clinical <- readRDS("../data/os_data.rds")

# read data
cluster_results <- readRDS("../results/euler_memberships.rds")
mutation_covariate_data <- readRDS("../data/aml_data.rds")

# merge features
clinical$group <- as.factor(cluster_results$clustermembership)
levels(clinical$group) <- LETTERS[1:9]

clinical$type <- mutation_covariate_data$Dx
clinical$gender <- mutation_covariate_data$Gender
clinical$age <- mutation_covariate_data$age



# Kaplan-Meier curve for groups
colourysdots <- c("#202020","#774411","#DDAA77","#ed2124","#114477","#CC99BB",
                           "#88CCAA","#117744","#77AADD")

os <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical)
km_plot <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = colourysdots, # personalized colours
  legend.labs = sort(as.character(unique(clinical$group))), 
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
); km_plot


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


# survival plot for each cancer type
clinical_aml <- as.data.frame(clinical)[clinical$type=="AML",]

os_aml <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical_aml)
ggsurvplot(
  os_aml,                     # survfit object with calculated statistics.
  data = clinical_aml,             # data used to fit survival curves.
  palette = colourysdots, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)

group_count <- clinical_aml %>% 
  group_by(group) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 5 times
filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_aml <- clinical_aml %>% 
  filter(group %in% filtered_groups)

os_aml <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = filtered_clinical_aml)
km_aml <- ggsurvplot(
  os_aml,                     # survfit object with calculated statistics.
  data = filtered_clinical_aml,             # data used to fit survival curves.
  palette = colourysdots[which(levels(clinical$group) %in% sort(as.character(unique(filtered_clinical_aml$group))))], # personalized colours
  legend.labs = sort(as.character(unique(filtered_clinical_aml$group))), 
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
); km_aml

# survival plot for each cancer type
clinical_mds <- as.data.frame(clinical)[clinical$type=="MDS",]

os_mds <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical_mds)
ggsurvplot(
  os_mds,                     # survfit object with calculated statistics.
  data = clinical_mds,             # data used to fit survival curves.
  palette = colourysdots, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

group_count <- clinical_mds %>% 
  group_by(group) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 5 times
filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_mds <- clinical_mds %>% 
  filter(group %in% filtered_groups)

os_mds <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = filtered_clinical_mds)
km_mds <- ggsurvplot(
  os_mds,                     # survfit object with calculated statistics.
  data = filtered_clinical_mds,             # data used to fit survival curves.
  palette = colourysdots[which(levels(clinical$group) %in% sort(as.character(unique(filtered_clinical_mds$group))))], # personalized colours
  legend.labs = sort(as.character(unique(filtered_clinical_mds$group))), 
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
); km_mds


# survival plot for each cancer type
clinical_cmml <- as.data.frame(clinical)[clinical$type=="CMML",]

os_cmml <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical_cmml)
ggsurvplot(
  os_cmml,                     # survfit object with calculated statistics.
  data = clinical_cmml,             # data used to fit survival curves.
  palette = colourysdots, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)

group_count <- clinical_cmml %>% 
  group_by(group) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 5 times
filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_cmml <- clinical_cmml %>% 
  filter(group %in% filtered_groups)

os_cmml <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = filtered_clinical_cmml)
km_cmml <- ggsurvplot(
  os_cmml,                     # survfit object with calculated statistics.
  data = filtered_clinical_cmml,             # data used to fit survival curves.
  palette = colourysdots[which(levels(clinical$group) %in% sort(as.character(unique(filtered_clinical_cmml$group))))], # personalized colours
  legend.labs = sort(as.character(unique(filtered_clinical_cmml$group))), 
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
); km_cmml

# survival plot for each cancer type
clinical_mpn <- as.data.frame(clinical)[clinical$type=="MPN",]

os_mpn <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical_mpn)
ggsurvplot(
  os_mpn,                     # survfit object with calculated statistics.
  data = clinical_mpn,             # data used to fit survival curves.
  palette = colourysdots, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)

group_count <- clinical_mpn %>% 
  group_by(group) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 5 times
filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_mpn <- clinical_mpn %>% 
  filter(group %in% filtered_groups)

os_mpn <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = filtered_clinical_mpn)
km_mpn <- ggsurvplot(
  os_mpn,                     # survfit object with calculated statistics.
  data = filtered_clinical_mpn,             # data used to fit survival curves.
  palette = colourysdots[which(levels(clinical$group) %in% sort(as.character(unique(filtered_clinical_mpn$group))))], # personalized colours
  legend.labs = sort(as.character(unique(filtered_clinical_mpn$group))), 
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
); km_mpn


# # ggarrange(km_aml$plot, km_aml_riskscores$plot, km_mds$plot, km_mds_riskscores$plot, widths = c(1,1), labels = c("a","b","c","d"))
# ggarrange(km_aml$plot, km_mds$plot, km_aml$table, km_mds$table, km_cmml$plot, 
#           km_mpn$plot, km_cmml$table, km_mpn$table, nrow = 4, ncol = 2, heights = c(1,0.5,1,0.5), widths = c(1,1), labels = c("a","b","","","c","d","",""))
# 
# 
# library("extrafont")
# loadfonts()
# pdf("~/Desktop/km_cancer_type.pdf", height = 10.7, width = 9,
#     family = "Arial", paper = "special", onefile = FALSE)
# # family = "Times New Roman", paper = "special", onefile = FALSE)
# op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
# ggarrange(km_aml$plot, km_mds$plot, km_aml$table, km_mds$table, km_cmml$plot, 
#           km_mpn$plot, km_cmml$table, km_mpn$table, nrow = 4, ncol = 2, heights = c(1,0.35,1,0.35), 
#           widths = c(1,1), labels = c("a","b","","","c","d","",""))
# par(op)
# dev.off()



most_frequent_group <- clinical_mpn %>%
  group_by(group) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n)) %>%
  slice_head(n = 1) %>%
  pull(group)

clinical_mpn_modified <- clinical_mpn %>%
  mutate(group = if_else(group == most_frequent_group, most_frequent_group, 'Other'))

os_mpn <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical_mpn_modified)
km_mpn_separate <- ggsurvplot(
  os_mpn,                     # survfit object with calculated statistics.
  data = clinical_mpn_modified,             # data used to fit survival curves.
  palette = c(colourysdots[which(levels(clinical$group) %in% sort(as.character(unique(filtered_clinical_mpn$group))))], "black"), # personalized colours
  legend.labs = sort(as.character(unique(clinical_mpn_modified$group))), 
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 1.7,
  pval.size =4,
  # title="My title",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
);km_mpn_separate

# ggarrange(km_aml$plot, km_aml_riskscores$plot, km_mds$plot, km_mds_riskscores$plot, widths = c(1,1), labels = c("a","b","c","d"))
ggarrange(km_aml$plot, km_mds$plot, km_aml$table, km_mds$table, km_cmml$plot,
          km_mpn_separate$plot, km_cmml$table, km_mpn_separate$table, nrow = 4, ncol = 2, heights = c(1,0.35,1,0.35), widths = c(1,1), labels = c("a","b","","","c","d","",""))

library("extrafont")
loadfonts()
pdf("~/Desktop/km_cancer_type.pdf", height = 10.7, width = 9,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
ggarrange(km_aml$plot, km_mds$plot, km_aml$table, km_mds$table, km_cmml$plot,
          km_mpn_separate$plot, km_cmml$table, km_mpn_separate$table, nrow = 4, ncol = 2, heights = c(1,0.35,1,0.35), 
          widths = c(1,1), labels = c("a","b","","","c","d","",""))

par(op)
dev.off()
