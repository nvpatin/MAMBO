library(ggplot2)
library("dplyr")
library("ggplot2")
library("mclust")
library("cluster")
library("tidyr")
library("randomForest")
# Example occurrence data -------------------------------------------------

num.samples <- 100
num.otus <- 500

num.reads <- matrix(
  sample(0:100, num.samples * num.otus, replace = TRUE),
  ncol = num.samples
)

# one group (1-50) has very high abundance of first 50 OTUs
num.reads[1:50, 1:50] <- num.reads[1:10, 1:10] + 1000
# one group (51-80) has medium high abundance of 51st to 60th OTUs
num.reads[51:60, 51:80] <- num.reads[51:60, 51:80] + 500

# Unsupervised Random Forest-------------------------------------------------

set.seed(131)                  #unsupervised RF
x.num.reads <- num.reads
unsup.num.reads <- randomForest(x.num.reads, ntree = 1000, proximity = TRUE)
unsup.mds.num.reads <- as.data.frame( cmdscale(1 - unsup.num.reads$proximity))

#Determine the number of clusters
Fing.mclust <- Mclust(unsup.mds.num.reads)
summary(Fing.mclust)
#plot(Fing.mclust) # plot results

#Set PAM cluster #
pamData <- pam(1 - unsup.num.reads$proximity, k = 3)
pamData$clustering
summary(Fing.mclust)
unsup.mds.num.reads$RFClassification <- Fing.mclust$classification
unsup.mds.num.reads$PAM <- pamData$clustering

#View(unsup.mds.num.reads)


PointRFClustersMclust <- ggplot(unsup.mds.num.reads, aes(x=V1, y=V2,  color = factor(RFClassification)) ) +  
  geom_point(size=4) + 
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
    panel.grid.major = element_line(colour = "white"),
    panel.background = element_rect(fill = "white", colour = "grey85")) +
  ggtitle("Random Forest Cluster") 

PointRFClustersPAM <- ggplot(unsup.mds.num.reads, aes(x=V1, y=V2,  color = factor(PAM)) ) +  
  geom_point(size=4) + 
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
    panel.grid.major = element_line(colour = "white"),
    panel.background = element_rect(fill = "white", colour = "grey85")) +
  ggtitle("Random Forest Cluster") 



