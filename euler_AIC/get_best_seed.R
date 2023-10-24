# get the best performing seed and corresponding DAGs + memberships

rm(list=ls())
setwd("/Users/frbayer/Documents/phd_main/projects/mda_analysis/euler_cluster")

assignprogressList <- list()
seeds <- c(1:79,81:100)
# seeds <- 1:100
for (i in seeds){
  assignprogressList[i] <- readRDS(paste0("euler_results/memberships_seed_",i,".rds"))$assignprogress
}

best_seed <- getBestSeed(assignprogressList)
best_res <- readRDS(paste0("euler_results/memberships_seed_",best_seed,".rds"))

saveRDS(best_res,"../results/euler_memberships.rds")
