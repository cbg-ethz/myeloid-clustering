# get the best performing seed and corresponding DAGs + memberships

rm(list=ls())
setwd("/Users/frbayer/Documents/phd_main/projects/mda_analysis/euler_cluster_number")
# setwd("/Users/frbayer/Documents/phd_main/projects/mda_analysis/euler_cluster_number_mutations")

library(ggplot2)
require(reshape2)
library(netClust)

cluster_range <- 5:30
chi_range <- 0:40*0.1

aic_matrix <- matrix(NA, nrow = length(cluster_range), ncol = length(chi_range))
kk <- 1
ii <- 1
for (k in cluster_range){
  ii <- 1
  for (i in chi_range){
    if(file.exists(paste0("euler_results/aic_k", k, "_chi", i, ".rds"))){
      aic_matrix[kk,ii] <- readRDS(paste0("euler_results/aic_k", k, "_chi", i, ".rds"))$testAIC
    }
    ii <- ii+1
  }
  kk <- kk+1
}



minK = 5; maxK = 30; chiVec = chi_range; AICrange = 200

# aic_res <- bestAICsearch_test(mut_cov_data, minK = minK, maxK = maxK, chiVec=chiVec, startseed = startseed, nIterations = nIterations, AICrange=AICrange)


### Prepare matrix for plotting
nrow <- maxK - minK + 1
ncol <- length(chiVec)
aics<-matrix(NA,nrow=nrow,ncol=ncol)

# aics <- c()
# for (hh in 1:length(aic_res)){
#   aics[hh] <- aic_res[[hh]]$testAIC
# }

aics <- aic_matrix

ks<-rep(0,ncol)
minaics<-rep(0,ncol)

minaics <- apply(aics, 2, min, na.rm=TRUE)

aics<-t(aics)
aics<-aics - minaics
topaics<-AICrange
aics[which(aics>topaics)]<-topaics
aicsscaled<-t(t(aics)/colSums(aics))
# heatmap(t(aics),Rowv=NA, Colv=NA,col = rainbow(256))

divergy<-aics
# divergy[which(is.na(divergy))]<-maxxy

rownames(divergy)<-chiVec
colnames(divergy)<-c(minK:maxK)

meltdivergy<-melt(divergy)

ggplot(data = meltdivergy, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()

middycol<-c(0.8,0.2,0)

ggheatmap<-ggplot(data = meltdivergy, aes(Var1, Var2, fill = value))+
  # ggheatmap<-ggplot(data = meltdivergy)+
  geom_tile() +
  xlab(expression(chi)) +
  ylab("k") +
  scale_fill_gradient2(high =rgb(0.98,0.98,1), low = rgb(0,0.35,0.8),
                       mid=rgb(0.49,0.665,0.9),space="Lab",na.value="grey75",
                       midpoint=topaics/2,limit = c(0,topaics), name="AIC\nchange\n") +
  scale_y_continuous(breaks=c(minK:maxK)) +
  theme_minimal() +
  theme(axis.title.x = element_text(vjust=-1),axis.title.y = element_text(angle=0,hjust=-0.5,vjust=0.505)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5,size = 20, hjust = 0.6),
        axis.text.y = element_text(angle = 0, vjust = 0.5,size = 20, hjust = 1),
        legend.text=element_text(size=20),
        axis.title=element_text(size=30),
        legend.title=element_text(size=24))+theme(legend.key.size = unit(2,"line")) +
  theme(plot.margin=unit(c(-0.3,-0.3,0.4,0.4),"cm"))

print(ggheatmap)

# pdf(paste("heatmapaic.pdf",sep=""), width=7.5, height=6, onefile=F, pointsize=10,  paper="special")

ggheatmap +
  theme(
    #axis.title.x = element_text("prune and reattach probability"),
    #axis.title.y = element_text("swap two nodes probability"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())







