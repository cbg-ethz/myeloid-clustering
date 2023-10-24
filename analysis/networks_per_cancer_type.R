# Learn networks for each cancer type

# Load required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggraph)
library(igraph)

# Clear the workspace
rm(list=ls())

# load data
cancer_type <- read.table("../data/categorical-mutationCovariate-matrix.txt")[,58]
mut_cov_data <- read.table("../data/categorical-mutationCovariate-matrix.txt")[,-58]
string_edgepmat <- as.matrix(read.table("../data/string-edgepmat.txt"))[-58,-58]
blacklist_edgepmat <- as.matrix(read.table("../data/blacklist-edgepmat.txt"))[-58,-58]

# legend for reference
# "AML" <- 0
# "CMML" <- 1
# "MDS" <- 2
# "MPN" <- 3

# # learn the cancer networks
# mut_cov_data$PB_AMC[6] <- 0 # because otherwise not present CMML
# mut_cov_data$BM_Gran[7] <- 0 # because otherwise not present in MPN
# maxfit <- list()
# bdepar=list(chi = 0.5, edgepf = 16)
# for (i in 1:4){
#   temp_mut_cov_data <- mut_cov_data[cancer_type==i-1,]
#   temp_mut_cov_data[1,colSums(temp_mut_cov_data)==0] <- 1
#   scorepar <- BiDAG::scoreparameters("bdecat", temp_mut_cov_data, edgepmat = string_edgepmat, bgnodes=c(58,59), bdepar = bdepar)
#   maxfit[[i]] <- BiDAG::iterativeMCMC(scorepar, verbose=FALSE, blacklist=blacklist_edgepmat)
# }
# saveRDS(maxfit, "../results/cancer_type_networks.rds")

maxfit <- readRDS("../results/cancer_type_networks.rds")

k_clust <- length(maxfit)

transform_DAG <- function(adj_matrix){
  # Transform the DAG to a CPDAG for a given variable selection
  adj_matrix_new <- adj_matrix
  temp_var_selection <- 1:46
  adj_matrix_new[temp_var_selection,temp_var_selection] <- as(pcalg::dag2cpdag(as(adj_matrix[temp_var_selection,temp_var_selection], "graphNEL")), "matrix")
  temp_var_selection <- 47:57
  adj_matrix_new[temp_var_selection,temp_var_selection] <- as(pcalg::dag2cpdag(as(adj_matrix[temp_var_selection,temp_var_selection], "graphNEL")), "matrix")
  
  return(adj_matrix_new)  
}

get_arrow_size <- function(mygraph, standard_size=2){
  # Identify directed and undirected edges
  edge_list <- as.data.frame(get.edgelist(mygraph))
  edge_list$undirected <- FALSE
  
  for (i in 1:nrow(edge_list)) {
    reversed_edge <- edge_list[i, c(2, 1)]
    if (any(edge_list$V1 == reversed_edge$V2 & edge_list$V2 == reversed_edge$V1)) {
      edge_list$undirected[i] <- TRUE
    }
  }
  
  # define arrow size
  arrow_size <- rep(standard_size, length(edge_list$undirected))
  arrow_size[edge_list$undirected] <- 0
  
  return(arrow_size)
}

# transform to CPDAGS where meaningful
for (hh in 1:4){
  maxfit[[hh]]$DAG <- transform_DAG(maxfit[[hh]]$DAG)
}

# plot with ggraph package
mycolor <- c("#e56d78", "#f7cac9",  "#708090", "#064273", "#92a8d1")

grp <- c(rep("Mutations",38), rep("Cytogenetic Alterations",8), rep("Blood Values",6), rep("Bone Marrow",5), rep("Adjusted Covariates",2))

# Compute by entropy
entropy <- function(cat.vect){
  # from https://stats.stackexchange.com/questions/221332/variance-of-a-distribution-of-multi-level-categorical-data
  px  = table(cat.vect)/length(cat.vect)
  lpx = log(px, base=2)
  ent = -sum(px*lpx)
  return(ent)
}

# top_entropy <- list()
# calculate overall entropy
total_entropies <- sapply(mut_cov_data, function(x) entropy(x))
# calculate relative entropy and select variables with the lowest entropy
n_net <- 16 # number of selected variables per cluster
all_var_selected <- c()
diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(mut_cov_data))
binary_frequency <- matrix(NA, nrow = k_clust, ncol = ncol(mut_cov_data))
for (nn in 1:k_clust){
  # selk_clustect data
  cluster_data <- mut_cov_data[cancer_type==nn-1,]
  # calculate entropy in clusters
  cluster_entropies <- sapply(cluster_data, function(x) entropy(x))
  # calculate relative entropy
  diff_entropies[nn,] <- cluster_entropies-total_entropies
  # calculate frequency
  binary_frequency[nn,] <- sapply(cluster_data, function(x) sum(x))
  # # select variables with lowest relative entropy
  # var_selection <- order(diff_entropies[nn,])[1:n_net]
  # all_var_selected <- unique(c(all_var_selected, var_selection))
  # select variables with highest frequency
  var_selection <- order(binary_frequency[nn,], decreasing = TRUE)[1:n_net]
  all_var_selected <- unique(c(all_var_selected, var_selection))
  # top_entropy[[nn]] <- var_selection
}

# set cols of genomic information and set frequ to zero elsewhere
binary_cols <- 1:46
binary_frequency[,-binary_cols] <- 0

# set node size
node_size_percentage <- 1-(diff_entropies[,]+abs(min(diff_entropies[,])))/(max(diff_entropies[,])+abs(min(diff_entropies[,])))
node_size <- 1.5+node_size_percentage*6

node_size_percentage_frequ <- (binary_frequency[,]+abs(min(binary_frequency[,])))/(max(binary_frequency[,])+abs(min(binary_frequency[,])))
node_size_frequ <- 1.5+node_size_percentage_frequ*6

node_size[,binary_cols] <- node_size_frequ[,binary_cols]

var_genes <- c(1:38)
var_cytogenetics <- c(39:46)
var_blood <- c(47:52)
var_bone <- c(53:57)
var_cov <- c(58:59)

all_var_selected <- all_var_selected[order(all_var_selected)]

# include age and sex
all_var_selected <- c(all_var_selected, 58, 59)

length(all_var_selected)

# include the most connected edges as well
k_clust <- 4
aa <- matrix()
adj_matrix <- maxfit[[1]]$DAG
aa <- colSums(as.matrix(adj_matrix)+t(as.matrix(adj_matrix)))
for (rr in 2:k_clust){
  adj_matrix <- maxfit[[rr]]$DAG
  aa <- aa+colSums(as.matrix(adj_matrix)+t(as.matrix(adj_matrix)))
}

n_connection <- 10

most_connected <- order(aa, decreasing = TRUE)[1:n_connection]
setdiff(most_connected, all_var_selected)

all_var_selected <- unique(c(most_connected, all_var_selected))
all_var_selected <- all_var_selected[order(all_var_selected)]

length(all_var_selected)

# specify node size of selected variables
node_size_selected <- node_size[,all_var_selected]

nn <- 1
n_included <- length(all_var_selected)
adj_matrix <- maxfit[[nn]]$DAG[all_var_selected,all_var_selected]
mygraph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="directed")

mycolor <- c("#708090", "#e56d78",  "#f7cac9", "#064273", "#92a8d1")

# grp <- c(rep("Mutations+Cytogenetics",n_included), rep("Blood Values",6), rep("Bone Marrow",5), rep("Conditioned Covariates",3))
grp <- c(rep("Mutations", length(intersect(var_genes, all_var_selected))), 
         rep("Cytogenetics", length(intersect(var_cytogenetics, all_var_selected))), 
         rep("Blood Values", length(intersect(var_blood, all_var_selected))), 
         rep("Bone Marrow", length(intersect(var_bone, all_var_selected))), 
         rep("Adjusted Covariates", length(intersect(var_cov, all_var_selected))))

# mycolor_vec <- c(rep("#e56d78",38), rep("#f7cac9",8), rep("#708090",6), rep("#064273",5), rep("#92a8d1",2))
mycolor_vec <- c(rep("#064273", length(intersect(var_genes, all_var_selected))), 
                 rep("#92a8d1", length(intersect(var_cytogenetics, all_var_selected))), 
                 rep("#e56d78", length(intersect(var_blood, all_var_selected))), 
                 rep("#f7cac9", length(intersect(var_bone, all_var_selected))), 
                 rep("#708090", length(intersect(var_cov, all_var_selected))))

mycolor_vec_alpha <- mycolor_vec
for (dd in 1:length(mycolor_vec)){
  mycolor_vec_alpha[dd] <- alpha(mycolor_vec[dd], node_size_percentage[nn,dd])
}

# Add label angle
number_of_bar=nrow(adj_matrix)
id = seq(1, nrow(adj_matrix))
angle= 360 * (id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
hjust <- ifelse(angle > 90 & angle<270, 1, 0)
angle <- ifelse(angle > 90 & angle<270, angle+180, angle)
name <- rownames(adj_matrix)
name[name=="cancer_type"] <- "Cancer Type"
name[name=="PB_ANC"] <- "ANC"
name[name=="PB_AMC"] <- "AMC"
name[name=="PB_WBC"] <- "WBC"
name[name=="PB_PLT"] <- "PLT"
name[name=="PB_HGB"] <- "HGB"
name[name=="PB_Blast"] <- "PB Blast"
name[name=="BM_Blast"] <- "BM Blast"
name[name=="BM_Mono"] <- "Mono"
name[name=="BM_Gran"] <- "Gran"
name[name=="BM_AuerRods"] <- "Auer Rod"
name[name=="BM_Dysplasias"] <- "Dysplasias"
Type <- as.factor(grp)

label_name <- name

p1 <- ggraph(mygraph, layout="circle") + 
  # geom_edge_link(edge_colour="black", edge_alpha=0.6, edge_width=0.4) +
  geom_edge_arc(arrow = arrow(length = unit(get_arrow_size(mygraph,2), 'mm')), 
                start_cap = circle(2.3, 'mm'),
                end_cap = circle(2, 'mm'), 
                edge_colour="black", edge_alpha=0.6, edge_width=0.4, aes(circular=TRUE)) +
  # geom_node_point(size=3.5, aes(color=Type), alpha=0.9) +
  geom_node_point(size=node_size_selected[nn,], aes(color=Type), alpha=0.9) +
  # scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=mycolor) +
  geom_node_text(aes(label=paste("    ",label_name,"    "), 
                     angle=angle, hjust=hjust), size=2.3, color="black") +
  theme_void() +
  theme(
    # legend.position="none",
    plot.margin=unit(c(0,0,0,0), "null"),
    panel.spacing=unit(c(0,0,0,0), "null")
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) + 
  coord_fixed(); p1



# plot all in a grid
p_list <- list()
k_clust <- length(maxfit)
for (jj in 1:k_clust){
  
  adj_matrix <- maxfit[[jj]]$DAG[all_var_selected,all_var_selected]
  mygraph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="directed")
  
  p_list[[jj]] <- ggraph(mygraph, layout="circle") +
    # geom_edge_link(edge_colour="black", edge_alpha=0.6, edge_width=0.4) +
    geom_edge_arc(arrow = arrow(length = unit(get_arrow_size(mygraph,1.5), 'mm'), type = 'closed'),
                  start_cap = circle(2.3, 'mm'),
                  end_cap = circle(2, 'mm'),
                  edge_colour="black", edge_alpha=0.6, edge_width=0.4, aes(circular=TRUE)) +
    # geom_node_point(size=3.5, aes(color=Type), alpha=0.9) +
    geom_node_point(size=node_size_selected[jj,], aes(color=Type), alpha=0.9) +
    # scale_size_continuous(range=c(0.5,8)) +
    scale_color_manual(values=mycolor) +
    geom_node_text(aes(label=paste("    ",label_name,"    "),
                       angle=angle, hjust=hjust), size=2.3, color="black") +
    theme_void() +
    theme(
      # legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5)) +
    coord_fixed() + theme(legend.text=element_text(size=13), legend.title = element_text(size=13)) 
  
}

# plot all on top of each other
ggarrange(plotlist=p_list, labels = c("AML", "CMML", "MDS", "MPN"),
          font.label = list(face = "plain"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

library("extrafont")
loadfonts()
pdf("~/Desktop/networks_cancer_type.pdf", height = 7, width = 10,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
ggarrange(plotlist=p_list, labels = c("AML", "CMML", "MDS", "MPN"),
          font.label = list(face = "plain"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")
par(op)
dev.off()

# plot all networks in one line

library("extrafont")
loadfonts()
pdf("~/Desktop/networks_cancer_type2.pdf", height = 3, width = 10,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
ggarrange(plotlist=p_list, labels = c("AML", "CMML", "MDS", "MPN"),
          font.label = list(face = "plain", size=13),
          ncol = 4, nrow = 1, common.legend = TRUE, legend="bottom")
par(op)
dev.off()

# plot the full networks

all_var_selected <- 1:59

# specify node size of selected variables
node_size_selected <- node_size[,all_var_selected]

nn <- 1
n_included <- length(all_var_selected)
adj_matrix <- maxfit[[nn]]$DAG[all_var_selected,all_var_selected]
mygraph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="directed")

mycolor <- c("#708090", "#e56d78",  "#f7cac9", "#064273", "#92a8d1")

# grp <- c(rep("Mutations+Cytogenetics",n_included), rep("Blood Values",6), rep("Bone Marrow",5), rep("Conditioned Covariates",3))
grp <- c(rep("Mutations", length(intersect(var_genes, all_var_selected))), 
         rep("Cytogenetics", length(intersect(var_cytogenetics, all_var_selected))), 
         rep("Blood Values", length(intersect(var_blood, all_var_selected))), 
         rep("Bone Marrow", length(intersect(var_bone, all_var_selected))), 
         rep("Adjusted Covariates", length(intersect(var_cov, all_var_selected))))

# mycolor_vec <- c(rep("#e56d78",38), rep("#f7cac9",8), rep("#708090",6), rep("#064273",5), rep("#92a8d1",2))
mycolor_vec <- c(rep("#064273", length(intersect(var_genes, all_var_selected))), 
                 rep("#92a8d1", length(intersect(var_cytogenetics, all_var_selected))), 
                 rep("#e56d78", length(intersect(var_blood, all_var_selected))), 
                 rep("#f7cac9", length(intersect(var_bone, all_var_selected))), 
                 rep("#708090", length(intersect(var_cov, all_var_selected))))

mycolor_vec_alpha <- mycolor_vec
for (dd in 1:length(mycolor_vec)){
  mycolor_vec_alpha[dd] <- alpha(mycolor_vec[dd], node_size_percentage[nn,dd])
}

# Add label angle
number_of_bar=nrow(adj_matrix)
id = seq(1, nrow(adj_matrix))
angle= 360 * (id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
hjust <- ifelse(angle > 90 & angle<270, 1, 0)
angle <- ifelse(angle > 90 & angle<270, angle+180, angle)
name <- rownames(adj_matrix)
name[name=="cancer_type"] <- "Cancer Type"
name[name=="PB_ANC"] <- "ANC"
name[name=="PB_AMC"] <- "AMC"
name[name=="PB_WBC"] <- "WBC"
name[name=="PB_PLT"] <- "PLT"
name[name=="PB_HGB"] <- "HGB"
name[name=="PB_Blast"] <- "PB Blast"
name[name=="BM_Blast"] <- "BM Blast"
name[name=="BM_Mono"] <- "Mono"
name[name=="BM_Gran"] <- "Gran"
name[name=="BM_AuerRods"] <- "Auer Rod"
name[name=="BM_Dysplasias"] <- "Dysplasias"
Type <- as.factor(grp)

label_name <- name

p1 <- ggraph(mygraph, layout="circle") + 
  # geom_edge_link(edge_colour="black", edge_alpha=0.6, edge_width=0.4) +
  geom_edge_arc(arrow = arrow(length = unit(get_arrow_size(mygraph,2), 'mm')), 
                start_cap = circle(2.3, 'mm'),
                end_cap = circle(2, 'mm'), 
                edge_colour="black", edge_alpha=0.6, edge_width=0.4, aes(circular=TRUE)) +
  # geom_node_point(size=3.5, aes(color=Type), alpha=0.9) +
  geom_node_point(size=node_size_selected[nn,], aes(color=Type), alpha=0.9) +
  # scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=mycolor) +
  geom_node_text(aes(label=paste("    ",label_name,"    "), 
                     angle=angle, hjust=hjust), size=2.3, color="black") +
  theme_void() +
  theme(
    # legend.position="none",
    plot.margin=unit(c(0,0,0,0), "null"),
    panel.spacing=unit(c(0,0,0,0), "null")
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) + 
  coord_fixed(); p1



# plot all in a grid
p_list <- list()
k_clust <- length(maxfit)
for (jj in 1:k_clust){
  
  adj_matrix <- maxfit[[jj]]$DAG[all_var_selected,all_var_selected]
  mygraph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="directed")
  
  p_list[[jj]] <- ggraph(mygraph, layout="circle") +
    # geom_edge_link(edge_colour="black", edge_alpha=0.6, edge_width=0.4) +
    geom_edge_arc(arrow = arrow(length = unit(get_arrow_size(mygraph,1.5), 'mm'), type = 'closed'),
                  start_cap = circle(2.3, 'mm'),
                  end_cap = circle(2, 'mm'),
                  edge_colour="black", edge_alpha=0.6, edge_width=0.4, aes(circular=TRUE)) +
    # geom_node_point(size=3.5, aes(color=Type), alpha=0.9) +
    geom_node_point(size=node_size_selected[jj,], aes(color=Type), alpha=0.9) +
    # scale_size_continuous(range=c(0.5,8)) +
    scale_color_manual(values=mycolor) +
    geom_node_text(aes(label=paste("    ",label_name,"    "),
                       angle=angle, hjust=hjust), size=2.3, color="black") +
    theme_void() +
    theme(
      # legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5)) +
    coord_fixed() + theme(legend.text=element_text(size=13), legend.title = element_text(size=13)) 
  
}

# plot all on top of each other
ggarrange(plotlist=p_list, labels = c("AML", "CMML", "MDS", "MPN"),
          font.label = list(face = "plain"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")


library("extrafont")
loadfonts()
pdf("~/Desktop/networks_ct_big_final.pdf", height = 9, width = 8,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
ggarrange(plotlist=p_list, labels = c("AML", "CMML", "MDS", "MPN"), 
          font.label = list(face = "plain"), 
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")
par(op)
dev.off()

