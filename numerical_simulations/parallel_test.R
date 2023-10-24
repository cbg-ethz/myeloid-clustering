# check speed-up

k_clust <- 4
n_vars <- 20
n_bg <- 10
n_samples <- NULL
bgedges <- "different"
equal_cpt_bg <- TRUE

set.seed(1)

fraction_list <- c()
difference_list <- c()

for (i in 1:10){

  # sample data
  sampled_results <- netClust:::sampleData(k_clust=k_clust, n_vars=n_vars, n_bg=n_bg, n_samples=n_samples,
                                           bgedges=bgedges, equal_cpt_bg=equal_cpt_bg)
  sampled_data <- sampled_results$sampled_data
  sampled_membership <- sampled_results$cluster_membership
  n_samples <- sampled_results$n_samples

  # sampled_results_list <- append(sampled_results_list, sampled_results)

  # clustering
  # correct_samples[ww,] <- netClust:::cluster_benchmark(sampled_data, sampled_membership, k_clust = k_clust,
  #                                                      n_bg = n_bg, n_vars = n_vars, n_rep = 1)


  ## cluster with covariate-adjusted framework
  start_time <- Sys.time()
  cluster_results1 <- netClust::get_clusters(sampled_data, k_clust = k_clust, n_bg = n_bg, EMseeds=i)
  end_time <- Sys.time()

  delta_time_parallel <- end_time-start_time


  ## cluster with covariate-adjusted framework
  start_time <- Sys.time()
  cluster_results1 <- netClust::netCluster(sampled_data, k_clust = k_clust, itLim = 50, n_bg = n_bg, EMseeds=i)
  end_time <- Sys.time()

  delta_time_nonparallel <- end_time-start_time


  fraction <- as.numeric(delta_time_nonparallel)/as.numeric(delta_time_parallel)

  difference <- as.numeric(delta_time_nonparallel)-as.numeric(delta_time_parallel)

  fraction_list <- c(fraction_list, fraction)
  difference_list <- c(difference_list, difference)

}



fraction_list1 <- fraction_list
difference_list1 <- difference_list

fraction 1.603422



