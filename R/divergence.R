# Calculate divergence between two probability density functions using the
# Kullback-Leibler distance measure.
# h - Returns matrix
# window - number of points to sample (defaults to anylength(h))
# count - number of observations to create
# filter - function to apply
# replace - whether to use replacement in bootstrapping
# divergence(sp500.subset, 25, filter=getCorFilter.RMT())
# divergence(sp500.subset, 25, filter=getCorFilter.Shrinkage())
# Can measure information (the default) or stability. Measuring stability will
# resample twice to get two forms of the correlation matrix.
divergence <- function(h, count, window=NULL, filter=getCorFilter.RMT(), 
  measure='information')
{
  fn <- paste('divergence', measure, sep='.')
  do.call(fn, list(h, count, window, filter))
}

divergence.information <- function(h, count, window, filter)
{
  if (is.null(window)) { window <- anylength(h) }
  # Convert to matrix to allow duplicates
  col.names <- colnames(h)
  row.names <- format(index(h), '%Y-%m-%d')
  h <- matrix(h, ncol=ncol(h))
  colnames(h) <- col.names
  rownames(h) <- row.names

  div <- function(junk, h.full)
  {
    h.window <- h.full[sample(index(h.full), window, replace=TRUE), ]
    c.sample <- cov2cor(cov.sample(h.window))
    c.model <- filter(h.window)

    divergence <- divergence.kl(c.sample, c.model)
    return(divergence)
  }
  ds <- sapply(1:count, div, h)

  theory <- divergenceLimit.kl(ncol(h), window)
  #cat("Theoretical divergence is",theory,"\n")

  return(c(mean=mean(ds, na.rm=TRUE), sd=sd(ds, na.rm=TRUE), limit=theory))
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

# The expected value of the divergence for random matrices (sample versus 
# true correlation matrix)
divergenceLimit.kl <- function(m, t=NULL)
{
  if (is.null(t))
  {
    t <- m[2]
    m <- m[1]
  }

  l <- t - m + 1
  0.5 * ( m * log(t/2) - sum(digamma((l:t)/2)) )
}

# plotDivergenceLimit.kl(100, 80:499, col='green', ylim=c(0,55))
# plotDivergenceLimit.kl(80, 80:499, col='orange', overlay=TRUE)
# plotDivergenceLimit.kl(40, 80:499, col='red', overlay=TRUE)
plotDivergenceLimit.kl <- function(m, t.range, ..., overlay=FALSE)
{
  ns <- rep(m,length(t.range))
  limit <- apply(matrix(c(ns, t.range), ncol=2), 1, divergenceLimit.kl)
  if (! overlay)
  {
    ylab <- 'Expected KL divergence'
    plot(t.range/m, limit, ylab=ylab, xlab='Q', type='l', ...)
  }
  else
  {
    lines(t.range/m, limit, type='l', ...)
  }

  invisible(limit)
}

# Limit for stability (distance between two sample correlation matrices)
stabilityLimit.kl <- function(m, t=NULL)
{ 
  if (is.null(t))
  {
    t <- m[2]
    m <- m[1]
  } 

  0.5 * m * (m+1) / (t - m - 1)
}

# Determine the stability of the filter.
divergence.stability <- function(h, count, window, filter)
{
  if (is.null(window)) { window <- anylength(h) }
  # Convert to matrix to allow duplicates
  h <- matrix(h, ncol=ncol(h))

  div <- function(junk, h.full)
  {
    h.window.1 <- h.full[sample(index(h.full), window, replace=TRUE), ]
    c.sample.1 <- cov2cor(cov.sample(h.window.1))
    h.window.2 <- h.full[sample(index(h.full), window, replace=TRUE), ]
    c.sample.2 <- cov2cor(cov.sample(h.window.2))

    divergence <- divergence.kl(c.sample.1, c.sample.2)
    return(divergence)
  }
  ds <- sapply(1:count, div, h)

  theory <- stabilityLimit.kl(ncol(h), window)
  #cat("Theoretical divergence is",theory,"\n")

  return(c(mean=mean(ds, na.rm=TRUE), sd=sd(ds, na.rm=TRUE), limit=theory))
}

