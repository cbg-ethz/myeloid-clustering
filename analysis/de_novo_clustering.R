# The following code describes the clustering of our data

# Network-based clustering of mutational and covariate data
library(clustNet)

# load mutationa and covariate data
mut_cov_data <- read.table("../data/binary-mutationCovariate-matrix.txt")[,-c(58:61)]
# the following data is optional, containing forbidden edges and priors from the STRING database
string_edgepmat <- as.matrix(read.table("../data/string-edgepmat-binary.txt"))[-c(58:61),-c(58:61)]
blacklist_edgepmat <- as.matrix(read.table("../data/blacklist-edgepmat-binary.txt"))[-c(58:61),-c(58:61)]

# define number of cluster-independent covariates (they need to be in the last cols of the matrix)
n_cov <- 2

# clustering
itLim <- 300
set.seed(330) # set seed

cluster_results <- clustNet::get_clusters(mut_cov_data, k_clust = 9, n_bg = n_cov,
                                quick = TRUE, EMseeds = 1,
                                bdepar=list(chi = 0.5, edgepf = 4),
                                edgepmat = string_edgepmat,
                                blacklist = blacklist_edgepmat)

# save the results
saveRDS(cluster_results, paste0("euler_results/memberships_seed_", seednumber, ".rds"))
