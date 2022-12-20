#! /usr/local/bin/Rscript 

rm(list = ls())
library(magrittr)
library(ggplot2)
source("/Users/nastassia.patin/GitHub/NOAA-NCAR-Hackathon/ranRelPct.R")

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


# Draw a sample -----------------------------------------------------------

# draw one sample of relative percentages and transpose matrix so rows are samples
ran1 <- t(ranRelPct(1, num.reads)[[1]])
# convert relative percentages to log-odds so euclidean distances make sense
lo.ran1 <- log(ran1 / (1 - ran1))
# compute pairwise euclidean distance of samples
lo.dist <- dist(lo.ran1)
# extract principal coordinates to data frame for plotting
lo.mds <- cmdscale(lo.dist) %>%
  as.data.frame() %>%
  setNames(c("PC1", "PC2"))


# Hierarchical clustering -------------------------------------------------

hclust.cl <- hclust(lo.dist)
lo.mds$hclust <- factor(cutree(hclust.cl, h = quantile(hclust.cl$height, 0.95)))

p1 <- ggdendro::ggdendrogram(hclust.cl)

p2 <- ggplot(lo.mds, aes(PC1, PC2)) +
  geom_point(aes(color = hclust)) +
  theme(legend.position = "none")

gridExtra::grid.arrange(p1, p2, nrow = 1, top = "Hierarchical clustering")


# K-means -----------------------------------------------------------------

kmeans.cl <- kmeans(lo.ran1, 3)
lo.mds$kmeans <- factor(kmeans.cl$cluster)
ggplot(lo.mds, aes(PC1, PC2)) +
  geom_point(aes(color = kmeans)) +
  ggtitle("K-Means") +
  theme(legend.position = "none")


# DBSCAN ------------------------------------------------------------------

minPts <- 4
knn <- data.frame(
  points = 1:num.samples,
  dist = sort(dbscan::kNNdist(lo.dist, k = minPts - 1))
)
eps <- 29.4

dbscan.cl <- dbscan::dbscan(lo.dist, eps = eps)
lo.mds$cluster <- factor(dbscan.cl$cluster)

p1 <- ggplot(knn, aes(points, dist)) +
  geom_line() +
  geom_hline(yintercept = eps, color = "red") +
  labs(
    x = paste0(minPts, "-NN distance"),
    y = "Points",
    title = paste("eps =", eps)
  )

p2 <- ggplot(lo.mds, aes(PC1, PC2)) +
  geom_point(aes(color = cluster)) +
  theme(legend.position = "none")

gridExtra::grid.arrange(p1, p2, nrow = 1, top = "DBSCAN")
