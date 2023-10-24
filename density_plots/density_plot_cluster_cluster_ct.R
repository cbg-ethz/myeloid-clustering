# Load required libraries
library(car)
library(ks)

# Clear environment and set working directory
rm(list=ls())

# Set parameters
excluded <- c(53:59) # focus on mutations + cytogenetics + blood values

# Load data
cancer_type <- read.table("../data/categorical-mutationCovariate-matrix.txt")[,58]
mut_cov_data <- read.table("../data/categorical-mutationCovariate-matrix.txt")[,-58][,-excluded]
string_edgepmat <- as.matrix(read.table("../data/string-edgepmat.txt"))[-58,-58][-excluded,-excluded]
blacklist_edgepmat <- as.matrix(read.table("../data/blacklist-edgepmat.txt"))[-58,-58][-excluded,-excluded]

cluster_results <- readRDS("../results/euler_memberships.rds")
k_clust <- length(cluster_results$DAGs)

# Calculate scores against clusters
scoresagainstclusters <- matrix(NA,dim(mut_cov_data)[1],k_clust)
for (k in 1:k_clust){
  allrelativeprobabs <- rep(0, dim(mut_cov_data)[1])
  allrelativeprobabs[cluster_results$clustermembership==k] <- 1
  
  scorepar <- BiDAG::scoreparameters("bdecat", mut_cov_data, weightvector = allrelativeprobabs, bdepar = bdepar)
  scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,cluster_results$DAGs[[k]][-excluded,-excluded],mut_cov_data, bdecatCvec = apply(mut_cov_data, 2, function(x) length(unique(x))))
}

# # Calculate scores against clusters
# scoresagainstclusters <- matrix(NA,dim(mut_cov_data)[1],k_clust)
# for (k in 1:k_clust){
#   temp_mut_cov_data <- mut_cov_data[cluster_results$clustermembership==k,]
#   scorepar <- BiDAG::scoreparameters("bdecat", temp_mut_cov_data, bdepar = bdepar)
#   scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,cluster_results$DAGs[[k]][-excluded,-excluded],mut_cov_data, bdecatCvec = apply(mut_cov_data, 2, function(x) length(unique(x))))
# }

# # there seem to be some NAs for the last cluster!!!
# scoresagainstclusters[is.na(scoresagainstclusters)] <- 2*min(scoresagainstclusters, na.rm = TRUE)

# Calculate divergence matrix
matsize <- nrow(mut_cov_data)
divergy <- matrix(0, matsize, matsize)
for(k1 in 1:(matsize-1)){
  P <- exp(scoresagainstclusters[k1,] - max(scoresagainstclusters[k1,]))
  P <- P/sum(P)
  for(k2 in (k1+1):matsize){
    Q <- exp(scoresagainstclusters[k2,] - max(scoresagainstclusters[k2,]))
    Q <- Q/sum(Q)
    M <- (P+Q)/2
    JS <- (P%*%log(P/M) + Q%*%log(Q/M))/2
    divergy[k1,k2] <- JS
    divergy[k2,k1] <- JS
  }
}

# Perform multidimensional scaling
reducey <- cmdscale(divergy, k=2)

plot(reducey, col=cancer_type+1)

# Load cluster and cancer data
# cluster_results <- readRDS("results/euler_memberships.rds")
# cluster_results$clustermembership <- cluster_results$clustermembership[[1]]
clustermembership <- as.factor(cluster_results$clustermembership)
levels(clustermembership) <- LETTERS[1:length(levels(clustermembership))]
levelly <- levels(clustermembership)

cancer_type <- as.factor(cancer_type)
levels(cancer_type) <- c("AML", "CMML", "MDS", "MPN")
levelly2 <- levels(cancer_type)

# Set colors and plotting parameters
colourys<-c("#202020","#774411","#DDAA77","#ed2124","#114477","#CC99BB",
                     "#88CCAA","#117744","#77AADD")
                     
colourys3<-alpha(colourys,0.3)
par(mar = c(0,0,0,0))

xlims<-c(1.08*min(reducey[,1]), max(reducey[,1]))
ylims<-c(min(reducey[,2]), max(reducey[,2]))

par(bty = 'n')
plot(1, type="n", axes=F, xlab="", ylab="",xlim=xlims,ylim=ylims, bty="n")

textheights<-(ylims[2]-ylims[1])*rev(1*((c(1:16)-1.5)/(16-1))) +ylims[1]
textxs<-xlims[1]

k_clust<-length(levelly)

for(k in 1:k_clust){
  selecteddots<-which(clustermembership==levelly[k])
  combdata<-reducey[selecteddots,]
  
  z <- kde(combdata)
  
  plot(z,display="filled.contour2",add=TRUE,cont=c(50),col=c("transparent",paste0(colourys[k],"66")), alpha=0.3, drawpoints=FALSE,drawlabels=FALSE,lwd=1.5, bty="n")  
  
  points(textxs,textheights[k],col=colourys[k],pch=19,cex=1.5)
  text(textxs,textheights[k],levelly[k],pos=4)
  
}


# select mutations
mutation_cols <- 1:46
mut_data <- mut_cov_data[,mutation_cols]

#TP53 check
tp53_index <- which(colnames(mut_data)=="TP53")
TP53only <- reducey[which((mut_data[,tp53_index]==1)&as.vector((rowSums(mut_data)==1)))[1],]
# text(TP53only[1], TP53only[2], "TP53 only")

#no muts check
nomuts <- reducey[which(as.vector(rowSums(mut_data)==0))[1],]

# plot edges for TP53 and no mutations
arrows(TP53only[1] + 0.05, TP53only[2] - 0.05, TP53only[1] + 0.005, TP53only[2] - 0.005, length=0.1)
text(TP53only[1] + 0.06, TP53only[2] - 0.085, paste0("TP53 only\n(",sum((mut_data[,1]==1)&as.vector((rowSums(mut_data)==1)))," samples)"))
arrows(nomuts[1] + 0.05, nomuts[2] + 0.05, nomuts[1] + 0.005, nomuts[2] + 0.005, length=0.1)
# text(nomuts[1] - 0.06, nomuts[2] - 0.085, paste0("No mutations\n(",sum(as.vector(rowSums(mut_data)==0))," samples)"))
text(nomuts[1] + 0.06, nomuts[2] + 0.085, paste0("No mutations and\n no cytogenetic abnormalities \n(",sum(as.vector(rowSums(mut_data)==0))," samples)"))


transparent_plot_plain <- recordPlot() 

colourys2 <- c("#DD7788", "#117777", "#771122", "#DDDD77")

# colourys2 <- alpha(RColorBrewer::brewer.pal(4, "Spectral"),0.3)

textxs<-xlims[2]*0.88

k_clust<-length(levelly2)

for(k in 1:k_clust){
  selecteddots<-which(cancer_type==levelly2[k])
  
  points(reducey[selecteddots,1],reducey[selecteddots,2],col=colourys2[k],pch=19,cex=0.5)
  
  # points(textxs,textheights[k],col=colourys[k],pch=19,cex=1.5)
  # text(textxs,textheights[k],levelly[k],pos=4)
  
}

transparent_plot <- recordPlot() 

pdf("~/Desktop/density_plot_cluster_cluster_ct.pdf", width = 9, height = 5, onefile=FALSE)
transparent_plot
dev.off()