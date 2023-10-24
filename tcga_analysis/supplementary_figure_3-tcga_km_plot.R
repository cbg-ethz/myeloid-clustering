# Prepare session, load packages
rm(list=ls())
library(survival)
library(RColorBrewer)
output <- read.table("data/tcga-clinical-information.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Load classification by mutation profile (cluster assignment)
cluster_res <- readRDS("results/adjust_memberships.rds")

output$cluster <- cluster_res$clustermembership

#rename some columns
colnames(output)[colnames(output)=="type"] <- "tissue"
colnames(output)[colnames(output)=="cluster"] <- "group"

# Change group from 1:22 to A:V
output$group <- as.factor(output$group)
levels(output$group) <- LETTERS[1:length(unique(output$group))]

# Merge clinical data with clustering
# clinical <- merge(output,cluster.membership, by="id", sort=FALSE)
clinical <- output
  
table(clinical$group)

clinical[,c('age','event','time')] <- lapply(clinical[,c('age','event','time')], as.numeric)
# Remove subtype from tissue
clinical[,c('stage','group','tissue')] <- lapply(clinical[,c('stage','group','tissue')], factor)

# Kaplan-Meier curve for 22 groups
# Colors A -- V

nb.cols <- length(unique(cluster_res$clustermembership))
# colourys <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
# set.seed(4)
# colourysdots <- sample(colourys,length(colourys))

colourysdots<-c("#202020","#771122","#AA4455","#DD7788","#774411","#AA7744",
                "#DDAA77","#777711","#AAAA44","#DDDD77","#117744","#44AA77",
                "#88CCAA","#117777","#44AAAA","#77CCCC","#114477","#4477AA",
                "#77AADD","#771155","#AA4488","#CC99BB")

par(mar = c(4.5,4.5,0.5,0.5))
plot(survfit(Surv(time = as.numeric(time)/365, event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', cex.lab=1.75, cex.axis=1.5, mark.time = T, xlim=c(0,24))
legend(x = 19, y = 1.01, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1.75, ncol = 2,title="Cluster", bg='white')

xs<-c(15,5.6,10,15,9.2,11.5,15.2,18.5,11.5,10.9,8.9,9.3,9.6,13.7,14.4,16.5,15,12.6,11.1,20.8,22.5,21.5)
ys<-c(0.86,0.27,0.34,0.645,0.385,0.445,0.115,0.125,0.505,0.03,0.425,0.2,0.545,0.1,0.28,0.17,0.4,0.03,0.32,0.28,0.125,0.41)
for(ii in 1:length(xs)){
  text(xs[ii],ys[ii],col=colourysdots[ii],LETTERS[ii],cex=1.5)
}

# pdf('~/Desktop/overall-survival-cluster-types-22.pdf', width=9, height=9)
# par(mar = c(4.5,4.5,0.5,0.5))
#     plot(survfit(Surv(time = as.numeric(time)/365, event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', cex.lab=1.75, cex.axis=1.5, mark.time = T, xlim=c(0,24))
#   legend(x = 19, y = 1.01, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1.75, ncol = 2,title="Cluster", bg='white')
#     
#     xs<-c(15,5.6,10,15,9.2,11.5,15.2,18.5,11.5,10.9,8.9,9.3,9.6,13.7,14.4,16.5,15,12.6,11.1,20.8,22.5,21.5)
#     ys<-c(0.86,0.27,0.34,0.645,0.385,0.445,0.115,0.125,0.505,0.03,0.425,0.2,0.545,0.1,0.28,0.17,0.4,0.03,0.32,0.28,0.125,0.41)
#     for(ii in 1:length(xs)){
#       text(xs[ii],ys[ii],col=colourysdots[ii],LETTERS[ii],cex=1.5)
#     }
#     
# dev.off()

library("extrafont")
loadfonts()
pdf('~/Desktop/overall-survival-cluster-types-22.pdf', width=15, height=8)
par(mar = c(4.5,4.5,0.5,0.5))
plot(survfit(Surv(time = as.numeric(time)/365, event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', cex.lab=1.75, cex.axis=1.5, mark.time = T, xlim=c(0,24))
legend(x = 21.5, y = 1.01, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1.75, ncol = 2,title="Cluster", bg='white')

xs<-c(15,5.6,10,15,9.2,11.5,15.2,18.5,11.5,10.9,8.9,9.3,9.6,13.7,14.4,16.5,15,12.6,11.1,20.8,22.5,21.5)
ys<-c(0.86,0.27,0.34,0.645,0.385,0.445,0.115,0.125,0.505,0.03,0.425,0.2,0.545,0.1,0.28,0.17,0.4,0.03,0.32,0.28,0.125,0.41)
for(ii in 1:length(xs)){
  text(xs[ii],ys[ii],col=colourysdots[ii],LETTERS[ii],cex=1.5)
}
dev.off()

