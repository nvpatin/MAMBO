rm(list = ls())

count <- 300
total <- 1000

a <- count + 1
b <- total - count + 1

curve(dbeta(x, a, b), 0, 1, 1000)
abline(v = count / total)


sampleRelPct <- function(x) {
  p.vec <- sapply(x, function(xi) {
    a <- xi + 1
    rbeta(1, a, sum(x) - a)
  })
  p.vec / sum(p.vec)
}

vec <- sample(0:5, 10, T)
mat <- t(replicate(1000, sampleRelPct(vec)))


vec / sum(vec)
apply(mat, 2, median)