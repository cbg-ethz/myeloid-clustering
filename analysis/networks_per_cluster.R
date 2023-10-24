# plot individual networks

# Load required libraries
library(ggplot2)
library(ggraph)
library(igraph)
library(ggpubr)
library(extrafont)
library(clustNet)

# Clear the workspace
rm(list=ls())

cluster_results <- readRDS("../results/euler_memberships.rds")
mut_cov_data <- read.table("../data/binary-mutationCovariate-matrix.txt")[,-c(58:61)]

# define useful functions
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
for (hh in 1:length(cluster_results$DAGs)){
  cluster_results$DAGs[[hh]] <- transform_DAG(cluster_results$DAGs[[hh]])
}

# plot with ggraph package
mycolor <- c("#ff727c", "#f7cac9",  "#708090", "#064273", "#92a8d1")
# mycolor <- c("#ff727c", "#f7cac9",  "#709090", "#064273", "#92a8d1")

grp <- c(rep("Mutations",38), rep("CG",8), rep("Blood Values",6), rep("Bone Marrow",5), rep("Adjusted Covariates",2))

# select cluster and corresponding DAG
nn <- 1 
adj_matrix <- cluster_results$DAGs[[nn]]
mygraph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="directed")

# Add label angle
number_of_bar=nrow(adj_matrix)
id = seq(1, nrow(adj_matrix))
angle= 360 * (id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
hjust <- ifelse(angle > 90 & angle<270, 1, 0)
angle <- ifelse(angle > 90 & angle<270, angle+180, angle)
name <- rownames(adj_matrix)
Type <- as.factor(grp)

k_clust <- length(cluster_results$DAGs)


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
n_net <- 15 # number of most frequent variables selected variables per cluster
all_var_selected <- c()
cluster_dim <- c()
diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(mut_cov_data))
binary_frequency <- matrix(NA, nrow = k_clust, ncol = ncol(mut_cov_data))
for (nn in 1:k_clust){
  # select data
  cluster_data <- mut_cov_data[cluster_results$clustermembership==nn,]
  cluster_dim[nn] <- dim(cluster_data)[1]
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

# normalize by cluster size
binary_frequency <- sweep(binary_frequency, 1, cluster_dim, FUN = "/")

binary_frequency_temp <- binary_frequency

# set cols of genomic information and set frequ to zero elsewhere
binary_cols <- 1:46
binary_frequency[,-binary_cols] <- 0

# set node size
node_size_percentage <- 1-(diff_entropies[,]+abs(min(diff_entropies[,])))/(max(diff_entropies[,])+abs(min(diff_entropies[,])))
node_size <- 1.5+node_size_percentage*6

node_size_percentage_frequ <- (binary_frequency[,]+abs(min(binary_frequency[,])))/(max(binary_frequency[,])+abs(min(binary_frequency[,])))
node_size_frequ <- 1.5+node_size_percentage_frequ*6

node_size[,binary_cols] <- node_size_frequ[,binary_cols]

# Find the maximum for each column and normalize per row for covariates 
max_col_values <- apply(binary_frequency_temp, 2, max)
node_size_percentage_frequ_temp <- sweep(binary_frequency_temp, 2, max_col_values, FUN = "/")
node_size_frequ_temp <- 1.5+node_size_percentage_frequ_temp*4

binary_cov_cols <- 47:58
node_size[,binary_cov_cols] <- node_size_frequ_temp[,binary_cov_cols]

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
aa <- matrix()
adj_matrix <- cluster_results$DAGs[[1]]
aa <- colSums(as.matrix(adj_matrix)+t(as.matrix(adj_matrix)))
for (rr in 2:k_clust){
  adj_matrix <- cluster_results$DAGs[[rr]]
  aa <- aa+colSums(as.matrix(adj_matrix)+t(as.matrix(adj_matrix)))
}

n_connection <- 10

most_connected <- order(aa, decreasing = TRUE)[1:n_connection]
setdiff(most_connected, all_var_selected)

all_var_selected <- unique(c(most_connected, all_var_selected))
all_var_selected <- all_var_selected[order(all_var_selected)]

length(all_var_selected)

# # remove cancer type
# all_var_selected <- all_var_selected[-which(all_var_selected==58)]

node_size_selected <- node_size[,all_var_selected]

nn <- 9
n_included <- length(all_var_selected)
adj_matrix <- cluster_results$DAGs[[nn]][all_var_selected,all_var_selected]
mygraph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="directed")

mycolor <- c("#708090", "#e56d78",  "#f7cac9", "#064273", "#92a8d1")
# mycolor <- c("#709090", "#e56d78",  "#f7cac9", "#064273", "#92a8d1")

grp <- c(rep("Mutations", length(intersect(var_genes, all_var_selected))), 
         rep("Cytogenetics", length(intersect(var_cytogenetics, all_var_selected))), 
         rep("Blood Values", length(intersect(var_blood, all_var_selected))), 
         rep("Bone Marrow", length(intersect(var_bone, all_var_selected))), 
         rep("Adjusted Covariates", length(intersect(var_cov, all_var_selected))))

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
name[name=="CG_5"] <- "-5/5q"
name[name=="CG_7"] <- "-7/7q"
name[name=="CG_8"] <- "+8/8q"
name[name=="CG_Y"] <- "-Y"
name[name=="CG_11"] <- "+11/11q"
name[name=="CG_inv16"] <- "inv16"
name[name=="CG_8_21"] <- "t(8;21)"
name[name=="CG_20"] <- "-20/20q"

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
k_clust <- length(cluster_results$DAGs)
for (jj in 1:k_clust){
  
  adj_matrix <- cluster_results$DAGs[[jj]][all_var_selected,all_var_selected]
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
ggarrange(plotlist=p_list, labels = paste("Cluster", LETTERS[1:k_clust]), 
          font.label = list(face = "plain"), 
          ncol = 3, nrow = 3, common.legend = TRUE, legend="bottom")

library("extrafont")
loadfonts()
pdf("~/Desktop/networks_final.pdf", height = 10, width = 10,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
ggarrange(plotlist=p_list, labels = paste("Cluster", LETTERS[1:k_clust]), 
          font.label = list(face = "plain"), 
          ncol = 3, nrow = 3, common.legend = TRUE, legend="bottom")
par(op)
dev.off()

# ------------

# plot full networks

all_var_selected <- 1:59

nn <- 9
n_included <- length(all_var_selected)
adj_matrix <- cluster_results$DAGs[[nn]][all_var_selected,all_var_selected]
mygraph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="directed")

# plot with ggraph package
grp <- c(rep("Mutations", length(intersect(var_genes, all_var_selected))), 
         rep("Cytogenetics", length(intersect(var_cytogenetics, all_var_selected))), 
         rep("Blood Values", length(intersect(var_blood, all_var_selected))), 
         rep("Bone Marrow", length(intersect(var_bone, all_var_selected))), 
         rep("Adjusted Covariates", length(intersect(var_cov, all_var_selected))))

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
name[name=="CG_5"] <- "-5/5q"
name[name=="CG_7"] <- "-7/7q"
name[name=="CG_8"] <- "+8/8q"
name[name=="CG_Y"] <- "-Y"
name[name=="CG_11"] <- "+11/11q"
name[name=="CG_inv16"] <- "inv16"
name[name=="CG_8_21"] <- "t(8;21)"
name[name=="CG_20"] <- "-20/20q"

Type <- as.factor(grp)

label_name <- name

p1 <- ggraph(mygraph, layout="circle") + 
  # geom_edge_link(edge_colour="black", edge_alpha=0.6, edge_width=0.4) +
  geom_edge_arc(arrow = arrow(length = unit(get_arrow_size(mygraph,2), 'mm')), 
                start_cap = circle(2.3, 'mm'),
                end_cap = circle(2, 'mm'), 
                edge_colour="black", edge_alpha=0.6, edge_width=0.4, aes(circular=TRUE)) +
  # geom_node_point(size=3.5, aes(color=Type), alpha=0.9) +
  geom_node_point(size=node_size[nn,], aes(color=Type), alpha=0.9) +
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
k_clust <- length(cluster_results$DAGs)
for (jj in 1:k_clust){
  
  adj_matrix <- cluster_results$DAGs[[jj]][all_var_selected,all_var_selected]
  mygraph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="directed")
  
  p_list[[jj]] <- ggraph(mygraph, layout="circle") +
    # geom_edge_link(edge_colour="black", edge_alpha=0.6, edge_width=0.4) +
    geom_edge_arc(arrow = arrow(length = unit(get_arrow_size(mygraph,1.5), 'mm'), type = 'closed'),
                  start_cap = circle(2.3, 'mm'),
                  end_cap = circle(2, 'mm'),
                  edge_colour="black", edge_alpha=0.6, edge_width=0.4, aes(circular=TRUE)) +
    # geom_node_point(size=3.5, aes(color=Type), alpha=0.9) +
    geom_node_point(size=node_size[jj,], aes(color=Type), alpha=0.9) +
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
    coord_fixed()+
    coord_fixed() + theme(legend.text=element_text(size=13), legend.title = element_text(size=13)) 
  
}

# plot all on top of each other
ggarrange(plotlist=p_list, labels = paste("Cluster", LETTERS[1:k_clust]), ncol = 3, nrow = 3, common.legend = TRUE, legend="bottom") #, font.label=list(color="black",size=9))

library("extrafont")
loadfonts()
pdf("~/Desktop/networks_big_final.pdf", height = 9.5, width = 8,
    family = "Arial", paper = "special", onefile = FALSE)
# family = "Times New Roman", paper = "special", onefile = FALSE)
op <- par(mar = c(5, 4, 0.05, 0.05) + 0.1)
ggarrange(plotlist=p_list, labels = paste("Cluster", LETTERS[1:k_clust]), 
          font.label = list(face = "plain"), 
          ncol = 3, nrow = 3, common.legend = TRUE, legend="bottom")
par(op)
dev.off()

