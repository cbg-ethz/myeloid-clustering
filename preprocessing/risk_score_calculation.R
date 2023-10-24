## calculate IPSSM for MDS ##

# load the ipssm library
library("ipssm")

# Clear the workspace
rm(list=ls())

# get example data
path.file <- system.file("extdata", "IPSSMexample.csv", package = "ipssm")
dd <- IPSSMread(path.file)

# read data
mutation_covariate_data <- readRDS("../data/aml_data.rds")

# rename columns according to ipssm library
names(mutation_covariate_data)[names(mutation_covariate_data) == 'UPI'] <- 'ID'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'PB_PLT'] <- 'PLT'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'PB_HGB'] <- 'HB'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'BM_Blast'] <- 'BM_BLAST'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'CG_5'] <- 'del5q'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'CG_7'] <- 'del7_7q'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'TP53'] <- 'TP53mut'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'MLL'] <- 'MLL_PTD'

# get intersect
var_intersect <- intersect(colnames(dd),colnames(mutation_covariate_data))

# define based on intersect
my_dd <- mutation_covariate_data[var_intersect]

# missing variables
var_missing <- setdiff(colnames(dd),var_intersect)
# missing: "complex"    ("CYTO_IPSSR") "del17_17p"  "TP53maxvaf" "TP53loh"    "GNB1"       "PPM1D"      "PRPF8"

# define NAs in missing variables
dd_temp <- dd[var_missing]
dd_temp[,] <- NA
dd_temp2 <- dd_temp[1:dim(my_dd)[1],]
rownames(dd_temp2) <- rownames(my_dd)

# get complete data
my_temp_final <- cbind(my_dd,dd_temp2[1:dim(my_dd)[1],])

# get score
my_dd.process <- IPSSMprocess(my_temp_final)
my_dd.res <- IPSSMmain(my_dd.process)
my_dd.annot <- IPSSMannotate(my_dd.res)

# store
mutation_covariate_data$IPSSM <- my_dd.annot$IPSSMcat_mean
# mutation_covariate_data$IPSSMscore_mean <- my_dd.annot$IPSSMscore_mean
mds_ids <- which(mutation_covariate_data$Dx=="MDS")
mutation_covariate_data$IPSSM[-mds_ids] <- NA

# save data
output_path <- "../data/aml_data_IPSSM_annotated.rds"
saveRDS(mutation_covariate_data, output_path)

## calculate ELN for AML ##
rm(list=ls())

# read data
mutation_covariate_data <- readRDS("../data/aml_data.rds")

colnames(mutation_covariate_data)


# rename columns according to ipssm library
names(mutation_covariate_data)[names(mutation_covariate_data) == 'CG_8_21'] <- 't(8;21)'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'CG_inv16'] <- 'inv(16)'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'CG_11'] <- 't(v;11)'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'CG_5'] <- '-5/del(5q)'
names(mutation_covariate_data)[names(mutation_covariate_data) == 'CG_7'] <- '-7/del(7q)'


var_missing <- setdiff(c("t(8;21)", "inv(16)", "NPM1", "CEBPA","FLT3", "t(9;11)","t(6;9)", "t(v;11)", "t(9;22)", "inv(3)", "-5/del(5q)",
                         "i(17)", "-7/del(7q)", "complex", "ASXL1", "BCOR",
                         "EZH2", "RUNX1", "SF3B1", "SRSF2", "STAG2", "U2AF1",
                         "ZRSR2", "TP53"),colnames(mutation_covariate_data))
# missing: t(9;11), t(6;9), t(9;22), inv(3), i(17), complex

# define function for calculation
calculate_risiko_ELN2022 <- function(dataframe) {
  
  favorable_genes <- c("t(8;21)", "inv(16)", "NPM1", "CEBPA")
  intermediate_genes <- c("FLT3", "t(9;11)")
  adverse_genes <- c("t(6;9)", "t(v;11)", "t(9;22)", "inv(3)", "-5/del(5q)",
                     "i(17)", "-7/del(7q)", "complex", "ASXL1", "BCOR",
                     "EZH2", "RUNX1", "SF3B1", "SRSF2", "STAG2", "U2AF1",
                     "ZRSR2", "TP53")
  
  for (i in 1:nrow(dataframe)) {
    # Skip if there is already a value in eln2022 for this row
    if (!is.na(dataframe$eln2022[i])) next
    
    gene_mutations <- colnames(dataframe)[which(dataframe[i,] == 1)]
    
    if (any(gene_mutations %in% favorable_genes)) {
      dataframe$eln2022[i] <- "Favorable (ELN2022)"
    } else if (any(gene_mutations %in% adverse_genes)) {
      dataframe$eln2022[i] <- "Adverse (ELN2022)"
    } else {
      dataframe$eln2022[i] <- "Intermediate (ELN2022)"
    }
  }
  
  return(dataframe)
}

mutation_covariate_data$eln2022 <- NA
mutation_covariate_data <- calculate_risiko_ELN2022(mutation_covariate_data)

# store
aml_ids <- which(mutation_covariate_data$Dx=="AML")
mutation_covariate_data$eln2022[-aml_ids] <- NA

# save data
output_path <- "../data/aml_data_eln2022_annotated.rds"
saveRDS(mutation_covariate_data, output_path)

## unify them in one file ##
eln_score <- readRDS("../data/aml_data_eln2022_annotated.rds")$eln2022
ipssm_score <- readRDS("../data/aml_data_IPSSM_annotated.rds")$IPSSM

# get rows of AML / MDS
aml_ids <- which(mutation_covariate_data$Dx=="AML")
mds_ids <- which(mutation_covariate_data$Dx=="MDS")

# label IPSSM better
ipssm_score <- as.factor(ipssm_score)
levels(ipssm_score) <- paste0(levels(ipssm_score)," (IPSSM)")
ipssm_score <- as.character(ipssm_score)

mutation_covariate_data <- readRDS("../data/aml_data_eln2022_annotated.rds")
mutation_covariate_data$ipssm <- ipssm_score
mutation_covariate_data$risk_scores <- NA
mutation_covariate_data$risk_scores[aml_ids] <- mutation_covariate_data$eln2022[aml_ids]
mutation_covariate_data$risk_scores[mds_ids] <- mutation_covariate_data$ipssm[mds_ids]

# save data
output_path <- "../data/risk_scores.rds"
saveRDS(mutation_covariate_data, output_path)


