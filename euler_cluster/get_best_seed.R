# get the best performing seed and corresponding DAGs + memberships
library(clustNet)

rm(list=ls())

assignprogressList <- list()
seeds <- c(1:500)
for (i in seeds){
  assignprogressList[[i]] <- readRDS(paste0("euler_results/memberships_seed_",i,".rds"))$assignprogress
}

best_seed <- clustNet:::getBestSeed(assignprogressList)
best_res <- readRDS(paste0("euler_results/memberships_seed_", best_seed,".rds"))

saveRDS(best_res,"../results/euler_memberships.rds")

