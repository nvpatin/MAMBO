rm(list = ls())

# Example occurrence data -------------------------------------------------

num.samples <- 5
num.otus <- 10

num.reads <- matrix(
  sample(0:30, num.samples * num.otus, replace = TRUE),
  ncol = num.samples
)

# Function to draw 'n' random relative percent matrices -------------------

ranRelPct <- function(n, num.reads) {
  coverage <- matrix(
    rep(colSums(num.reads), each = nrow(num.reads)),
    ncol = ncol(num.reads)
  )

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

# This is the mode (expected value) of the distribution
rel.pct <- t(t(num.reads) / colSums(num.reads))

# One random draw
ran.1 <- ranRelPct(1, num.reads)

# Ten random draws
ran.10 <- ranRelPct(10, num.reads)
