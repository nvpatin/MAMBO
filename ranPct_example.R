rm(list = ls())

# Example occurrence data -------------------------------------------------

num.samples <- 5
num.otus <- 10
coverage.range <- c(10, 1000)

# matrix of coverage
coverage <- matrix(
  sample(
    min(coverage.range):max(coverage.range),
    num.samples * num.otus,
    replace = TRUE
  ),
  ncol = num.samples
)

# matrix of reads (limited by coverage
num.reads <- matrix(
  sapply(coverage, function(x) sample(0:x, 1)),
  ncol = num.samples
)


# Function to draw 'n' random % occurrence matrices -----------------------

ranPctReads <- function(n, num.reads, coverage) {
  # fit beta shape parameters
  beta.params <- array(
    c(num.reads + 1, coverage - num.reads + 1),
    dim = c(dim(num.reads), 2)
  )

  # draw n random matrices from beta distribution
  replicate(n, {
    apply(beta.params, c(1, 2), function(x) rbeta(1, x[1], x[2]))
  }, simplify = FALSE)
}


# Test it -----------------------------------------------------------------

occ.mode <- num.reads / coverage
ran.1 <- ranPctReads(1, num.reads, coverage)
ran.10 <- ranPctReads(10, num.reads, coverage)
