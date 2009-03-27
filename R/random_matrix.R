# Look at rmetrics.org
#library('PerformanceAnalytics')
library(zoo)
library(quantmod)

# This is the default filter
getCorFilter.RMT <- function(breaks=NULL, hint=c(4,1))
{
  # h is a zoo object and will be t x m
  #function(h) { filter.RMT(t(h), breaks, hint) }
  function(h) return(filter.RMT(h, breaks, hint))
}

# This acts as a control case with no cleaning
getCorFilter.raw <- function()
{
  function(h) return(cor.empirical(h))
}

# Optimize a portfolio to minimize risk using RMT
# Params
#  h: a matrix of returns. Normalization occurs within the process
#  window: window size for calculating optimization. There would then be
#    t - window weights vectors returned.
#  breaks: MP density histogram breaks
#  hint: Initial parameters for MP fit optimization
# Returns
#  A weights vector or a matrix of weights over t - window steps
# Usage:
#  h <- some.returns.matrix # Dims: m x t; not normalized
#  ws <- optimizePortfolio.RMT(h)
# Example:
#  c.weights.260 <- optimizePortfolio.RMT(spx.small.history, window=260)
#  c.stat.260 <- plotPerformance(spx.small.history, c.weights.260)
#  # Compare to equal weighted portfolio
#  s.weights <- rep(1/50, 50)
#  s.stat <- plotPerformance(spx.small.history[,260:ncol(spx.small.history)], s.weights, TRUE, color="grey")
# NOTE: This is superceded by optimizePortfolio
optimizePortfolio.RMT <- function(h, window=NULL, breaks=NULL, hint=c(4,1))
{
  if (! 'zoo' %in% class(h))
  { cat("WARNING: Zoo objects are preferred for dimensional safety\n") }

  if (is.null(window)) { window <- nrow(h) }

  my.rmt <- function(h)
  {
    # Internally we work on m x t
    c.denoised <- filter.RMT(t(h), breaks, hint)

    if (logLevel() > 0) { cat("Optimizing portfolio weights\n") }
    p.optimize(h, c.denoised)
  }

  rollapply(h, window, my.rmt, by.column=FALSE, align='right')
}

# Transition in progress to TxM - filter.RMT now takes TxM xts object
filter.RMT <- function(h, breaks=NULL, hint=c(4,1))
{
  log.level <- logLevel()

  if (log.level > 1) { cat("Calculating eigenvalue distribution\n") }
  #mp.hist <- mp.density.hist(h, breaks=breaks, showHist=usePlots())
  mp.hist <- mp.density.kernel(h)

  mp.params <- optim(hint, mp.fit.kernel(mp.hist))
  mp.Q <- mp.params$par[1]
  mp.sigma <- mp.params$par[2]
  if (log.level > 0) { cat("Fit to Marcenko-Pastur with Q",mp.Q,"and sigma",mp.sigma,"\n") }


  if (usePlots())
  {
    if (anylength(dev.list()) > 0)
      old.par <- par(new=TRUE, ann=FALSE, yaxt='n', col='#44549C', bty='n')
    else
      old.par <- par(ann=FALSE, yaxt='n', col='#44549C', bty='n')
    rho <- mp.theory(mp.Q, mp.sigma)
    par(old.par)
  }

  if (log.level > 1) { cat("Cleaning eigenvalues\n") }
  # Use sigma^2 = 1 - lambda.1 / n
  #lambda.plus <- mp.eigen.max(mp.Q, mp.sigma)
  lambda.1 <- mp.hist$values[1]
  sigma.2 <- sqrt(1 - lambda.1/length(mp.hist$values))
  lambda.plus <- mp.eigen.max(mp.Q, sigma.2)
  if (log.level > 0)
  {
    cat("Upper cutoff (lambda.max) is",lambda.plus,"\n")
  }
  if (log.level > 1)
  {
    cat("Greatest eigenvalue is",lambda.1,"\n")
    cat("Variance is",sigma.2,"\n")
  }
  e.clean <- clean.bouchaud(t(h), mp.hist$values, lambda.plus)

  if (log.level > 1) { cat("Cleaning correlation matrix\n") }
  denoise(e.clean, mp.hist$vectors, t(h))
}

# h TxM zoo returns matrix
cor.empirical <- function(h)
{
  # Normalize returns
  ns <- r.normalize(h)

  # Calculate the correlation matrix
  # E = H H'
  t <- nrow(ns)
  e <- t(ns) %*% ns / t

}

##------------------------ MARCENKO-PASTUR FUNCTIONS ------------------------##
# Calculate the Marcenko-Pastur density
# h is now TxM
# t should be selected based on t = Q * m, where Q defines the shape of
# the Marcenko-Pastur density
# Example
#   h <- array(rnorm(t*m), c(t,m))
#   ee <- mp.density(h)[0]
#   ee <- mp.density(h, breaks=seq(0.01,3.01,0.02))
# Returns hist with eigenvalues attached
mp.density.hist <- function(h, showHist=TRUE, breaks=NULL, cutoff=0.01)
{
  e <- cor.empirical(h)

  # Calculate eigenvalues
  lambda <- eigen(e, symmetric=TRUE, only.values=FALSE)
  xbounds = c(0, floor(lambda$values[1] + 2))
  #cat("Eigenvalue bounds:",xbounds,"\n")
  if (is.null(breaks))
  {
    # Empirical value that gives okay breaks
    step = max(0.1, 20 / nrow(h)^1.5)
    breaks = seq(0, (xbounds[2]), step)
    #cat("Max:",xbounds[2],"\nRows:",step,"\n","Breaks:",breaks,"\n")
  }
  #cat("Values:",lambda$values)
  label <- paste('Eigenvalues (step size=',step,')',sep='')
  par(bg='#DDDDDD', cex.axis=0.6, cex.lab=0.8, col='#9599BB')

  hist <- hist(lambda$values[lambda$values > cutoff], 
    breaks=breaks, main=NA, xlab=label, freq=FALSE, xlim=c(0,6), plot=showHist)

  # If properly normalized, these should sum to m
  hist$values <- lambda$values
  hist$vectors <- lambda$vectors
  hist
}

mp.density.kernel <- function(h, ...)
{
  e <- cor.empirical(h)

  # Calculate eigenvalues
  lambda <- eigen(e, symmetric=TRUE, only.values=FALSE)
  ds <- density(lambda$values, ...)
  ds$values <- lambda$values
  ds$vectors <- lambda$vectors
  plot(ds, xlim=c(0,6))
  return(ds)
}

# Calculate and plot the theoretical density distribution
mp.theory <- function(Q, sigma, lambda=NULL, steps=200)
{
  # Plot a range of values
  if (is.null(lambda)) { lambda <- mp.lambdas(Q,sigma,steps) }
  rho <- mp.rho(Q,sigma, lambda)

  if (length(lambda) > 1)
  {
    l.min <- mp.eigen.min(Q,sigma)
    l.max <- mp.eigen.max(Q,sigma)
    xs <- seq(round(l.min-1), round(l.max+1), (l.max-l.min)/steps)
    plot(xs,rho, xlim=c(0,6), type='l')
  }
  rho
}


# Generate eigenvalues for theoretical Marcenko-Pastur distribution
mp.lambdas <- function(Q,sigma, steps)
{
  l.min <- mp.eigen.min(Q,sigma)
  l.max <- mp.eigen.max(Q,sigma)

  log.level <- logLevel()
  if (log.level > 2)
  {
    cat("min eigenvalue:",l.min,"\n")
    cat("max eigenvalue:",l.max,"\n")
  }

  evs <- seq(round(l.min-1), round(l.max+1), (l.max-l.min)/steps)
  evs[evs < l.min] <- l.min
  evs[evs > l.max] <- l.max

  if (log.level > 5)
  {
    #cat("x labels: ",xs,"\n")
    cat("eigenvalues: ",evs,"\n")
  }
  evs
}

# This provides the density of the eigenvalues
# lambda can be a vector of eigen values or a single eigen value
mp.rho <- function(Q,sigma, lambda)
{
  l.min <- mp.eigen.min(Q,sigma)
  l.max <- mp.eigen.max(Q,sigma)

  k <- (Q / 2*pi*sigma^2)
  rho <- k * sqrt(pmax(0, (l.max-lambda)*(lambda-l.min)) ) / lambda
  rho
}

# Get maximum eigenvalue specified by Marcenko-Pastur
mp.eigen.max <- function(Q,sigma) { sigma^2 * (1 + sqrt(1/Q))^2 }

# Get minimum eigenvalue specified by Marcenko-Pastur
mp.eigen.min <- function(Q,sigma) { sigma^2 * (1 - sqrt(1/Q))^2 }

# Fit the appropriate MP curve to the data. This will estimate Q and sigma.
# hist is a histogram object
# To simulate a fit, generate a random returns matrix:
#   h <- getRandomMatrix(20, 100)
#   hist <- mp.density(h)
#   hint <- c(6,2)
#   o <- optim(hint, mp.fit(hist))
#
#   load('RandomMatrix/spxReturnsQ5.rData')
#   spx.small <- t(as.matrix(spx.retq5[1:250,][,1:50]))
# NOTE: The scale factor needs to be exposed to adjust the fit, or just use
# a different function.
# Also, the cutoff in the cleaning can be adjusted to improve the filtering
mp.fit.hist <- function(hist)
{
  really.big <- 10000000000000
  log.level <- logLevel()
  fn <- function(ps)
  {
    Q <- ps[1]
    sigma <- ps[2]
    if (Q <= 0)
    { 
      if (log.level > 1) { cat("Kicking out Q<0:",Q,"sigma:",sigma,"\n") }
      return(really.big)
    }
    # Empirically, sigmas below 0.6 are unrealistic
    if (sigma <= 0.6)
    {
      if (log.level > 1) { cat("Kicking out sigma<0.6:",Q,"sigma:",sigma,"\n") }
      return(really.big)
    }

    l.plus <- mp.eigen.max(Q,sigma)
    rhos <- mp.rho(Q,sigma, hist$mids)

    # Just use some very large number to prevent it from being used as optimal
    # score
    if (max(rhos) == 0)
    {
      if (log.level > 1) { cat("Kicking out max(rho)==0:",Q,"sigma:",sigma,"\n") }
      return(really.big)
    }

    # Normalize based on amount of density below MP upper limit
    norm.factor <- sum(hist$density[hist$mids <= l.plus]) * 0.1
    if (log.level > 1) { cat("Using Q:",Q,"sigma:",sigma,"l+:",l.plus,"\n") }

    # Scale distance to be inline with calculated densities and histogram
    # This is a bit of hand-waving to get the best fit
    #scale <- max(rhos) / max(hist$density) + 1
    scale <- max(rhos) / max(hist$density) + 0.25
    if (log.level > 5) { cat("rhos:",rhos,"\n") }
    if (log.level > 2) { cat("scale:",scale,"\n") }

    dy <- (rhos - (hist$density * scale)) / norm.factor
    dist <- as.numeric(dy %*% dy)
    if (log.level > 5) { cat("dy:",dy,"\n\n") }
    if (log.level > 2) { cat("dist:",dist,"\n\n") }

    dist
  }
  fn
}

mp.fit.kernel <- function(hist)
{
  really.big <- 10000000000000
  log.level <- logLevel()
  fn <- function(ps)
  {
    Q <- ps[1]
    sigma <- ps[2]
    if (Q <= 0)
    { 
      if (log.level > 1) { cat("Kicking out Q<0:",Q,"sigma:",sigma,"\n") }
      return(really.big)
    }
    # Empirically, sigmas below 0.6 are unrealistic
    if (sigma <= 0.6)
    {
      if (log.level > 1) { cat("Kicking out sigma<0.6:",Q,"sigma:",sigma,"\n") }
      return(really.big)
    }

    l.plus <- mp.eigen.max(Q,sigma)
    rhos <- mp.rho(Q,sigma, hist$x)

    # Just use some very large number to prevent it from being used as optimal
    # score
    if (max(rhos) == 0)
    {
      if (log.level > 1) { cat("Kicking out max(rho)==0:",Q,"sigma:",sigma,"\n") }
      return(really.big)
    }

    # Normalize based on amount of density below MP upper limit
    norm.factor <- sum(hist$y[hist$x <= l.plus]) * 0.1
    if (log.level > 1) { cat("Using Q:",Q,"sigma:",sigma,"l+:",l.plus,"\n") }

    # Scale distance to be inline with calculated densities and histogram
    # This is a bit of hand-waving to get the best fit
    #scale <- max(rhos) / max(hist$density) + 1
    scale <- max(rhos) / max(hist$y) + 0.25
    if (log.level > 5) { cat("rhos:",rhos,"\n") }
    if (log.level > 2) { cat("scale:",scale,"\n") }

    dy <- (rhos - (hist$y * scale)) / norm.factor
    dist <- as.numeric(dy %*% dy)
    if (log.level > 5) { cat("dy:",dy,"\n\n") }
    if (log.level > 2) { cat("dist:",dist,"\n\n") }

    dist
  }
  fn
}

##---------------------------- RETURNS FUNCTIONS ----------------------------##
# Create a random returns matrix with population m (num securities) and 
# observations t (num returns per security)
# TODO: Make this zoo compatible (basically create a dummy index)
getRandomMatrix <- function(m, t)
{
  array(rnorm(m*t), c(t,m))
}

# Normalizes a returns matrix such that Var[xit] = 1.
# Assumes TxM, m population, t observations
r.normalize <- function(h) { apply(h, 2, function(x) x / sd(x)) }


##------------------------ CORRELATION MATRIX FUNCTIONS ---------------------##
# Normalize eigenvalues below some threshold lambda
# This is based on the Bouchaud cleaning trick
# Params:
#  h: m x t returns matrix
#  lambdas: vector of eigenvalues
#  lambda.plus: The upper bound that is used to cutoff noise. The default value
#    is a generic upper lambda bound of the MP distribution.
# Returns:
#  cleaned eigenvalues with noise values normalized
clean.bouchaud <- function(h, lambdas, lambda.plus=1.6)
{
  #m <- nrow(h)
  ##lam.p <- lambda.plus * 1.1
  #lam.p <- lambda.plus + 0.3
  #lam.m <- lambda.plus * 0.9
  #lambda.star <- sum(lambdas[lambdas<lam.p]) / length(lambdas[lambdas<lam.p])
  #clean = numeric(m)
  #for (i in 1:m)
  #{
  #  if (lambdas[i] > lam.m) { clean[i] <- lambdas[i] }
  #  else { clean[i] <- lambda.star }
  #}
  #clean
  avg <- mean(lambdas[lambdas < lambda.plus])
  lambdas[lambdas < lambda.plus] <- avg
  return(lambdas)
}

# Undo diagonalization of sample correlation matrix
# This takes flattened eigenvalues and returns a new cleaned correlation matrix
# Params:
#  e.values: Cleaned eigenvalues
#  e.vectors: Eigenvectors of correlation matrix of normalized returns
#  h: non-normalized returns matrix (only used for labels)
denoise <- function(e.values, e.vectors, h=NULL)
{
  #c.clean <- e.vectors %*% diag(e.values) %*% t(e.vectors)
  #diag(c.clean) <- 1
  ## These cleaned eigen values aren't really used
  ##c.denoised <- eigen(c.clean, symmetric=TRUE, only.values=FALSE)
  #if (! is.null(h))
  #{
  #  rownames(c.clean) <- rownames(h)
  #  colnames(c.clean) <- rownames(h)
  #}
  #c.clean
  c.clean <- e.vectors %*% diag(e.values) %*% t(e.vectors)
  diags <- diag(c.clean) %o% rep(1, nrow(c.clean))
  return(c.clean / sqrt(diags * t(diags)))
}


