# Calculate divergence between two probability density functions using the
# Kullback-Leibler distance measure.
divergence <- function(returns, size, count, fun, ...)
{
  # Resample count times
  bootstraps <- list()
  for (i in 1:count)
  {
    samples <- sample(1:nrow(returns), size)
    bootstraps[[i]] <- returns[samples,]
  }

  fn <- function(h)
  {
    cor.sample <- cov2cor(cov.sample(h))
    cor.filter <- fun(h)
    # Later make this dynamic in case someone comes up with something else
    divergence.kl(cor.sample, cor.filter, ...)
  }
  lapply(bootstraps, fn)
}

# Measuring information compares sample correlation matrix with filtered
# correlation matrix
# Measuring stability averages all permutations of the KL divergence of two
# instances of the filtered correlation matrix
divergence.kl <- function(sigma.1, sigma.2)
{
  term.1 <- log(det(sigma.2) / det(sigma.1))
  term.2 <- sum(diag(solve(sigma.2) %*% sigma.1))
  0.5 * (term.1 + term.2 - nrow(sigma.1))
}

# The expected value of the divergence for random matrices
divergenceLimit.kl <- function(n, T=NULL)
{
  if (is.null(T))
  {
    T <- n[2]
    n <- n[1]
  }

  l <- T - n + 1
  0.5 * ( n * log(T/2) - sum(digamma((l:T)/2)) )
}

# plotDivergenceLimit.kl(100, 80:499, col='green', ylim=c(0,55))
# plotDivergenceLimit.kl(80, 80:499, col='orange', overlay=TRUE)
# plotDivergenceLimit.kl(40, 80:499, col='red', overlay=TRUE)
plotDivergenceLimit.kl <- function(n, T.range, ..., overlay=FALSE)
{
  ns <- rep(n,length(T.range))
  limit <- apply(matrix(c(ns, T.range), ncol=2), 1, divergenceLimit.kl)
  if (! overlay)
  {
    ylab <- 'Expected KL divergence'
    plot(T.range/n, limit, ylab=ylab, xlab='Q', type='l', ...)
  }
  else
  {
    lines(T.range/n, limit, type='l', ...)
  }

  invisible(limit)
}

