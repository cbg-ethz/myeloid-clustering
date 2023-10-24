# Data processing string to edgepmat
rm(list=ls())

library("igraph")

# load mutation and covariate matrix
mutCovData <- read.table("../data/binary-mutationCovariate-matrix.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# load STRING graph and remove irrelevant variables
string <- as.matrix(read.table("../data/string-network.txt"))
irrelevantV <- setdiff(colnames(string),colnames(mutCovData))
stringIG <- igraph::graph_from_adjacency_matrix(string)
stringIG <- delete_vertices(stringIG, irrelevantV)

# create empty graph including all variabels
Nvar <- length(colnames(mutCovData))
edgepmatAll <- matrix(0,Nvar,Nvar)
colnames(edgepmatAll) <- colnames(mutCovData)
rownames(edgepmatAll) <- colnames(mutCovData)
edgepmatAllIG <- igraph::graph_from_adjacency_matrix(edgepmatAll)

# combine graphs (empty and STRING)
combinedIG <- igraph::union(edgepmatAllIG,stringIG)
edgepmat <- as.matrix(as_adjacency_matrix(combinedIG))
edgepmat2 <- 2-edgepmat
# check whether sorting is correct
colnames(edgepmat2)==colnames(mutCovData)

# save edgepmat
write.table(edgepmat2, "../data/string-edgepmat-binary.txt")

# create blacklist
edgepmat2_alternative <- edgepmat2
edgepmat2_alternative[,] <- 0
edgepmat2_alternative[47:61,] <- 1
edgepmat2_alternative[,47:57] <- 0
write.table(edgepmat2_alternative, "../data/blacklist-edgepmat-binary.txt")

