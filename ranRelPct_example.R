rm(list = ls())
source("ranRelPct.R")

# Example occurrence data -------------------------------------------------

num.samples <- 5
num.otus <- 10

num.reads <- matrix(
  sample(0:30, num.samples * num.otus, replace = TRUE),
  ncol = num.samples
)


# Test it -----------------------------------------------------------------

# This is the mode (expected value) of the distribution
rel.pct <- t(t(num.reads) / colSums(num.reads))

# One random draw
ran.1 <- ranRelPct(1, num.reads)

# Ten random draws
ran.10 <- ranRelPct(10, num.reads)


