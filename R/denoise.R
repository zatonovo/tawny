# RandomMatrixFilters consist of three components:
# . cor.fn - A function
# . cutoff.fn - A function
# . clean.fn - A function
create.RandomMatrixFilter <- function(T, cor.fn=cor.empirical,
  cutoff.fn=cutoff, clean.fn=cor.clean, ...)
{
  list(cor.fn=cor.fn, cutoff.fn=cutoff.fn, clean.fn=clean.fn, ...)
}

# Specify market as a symbol if you want to shrink on residuals only
create.ShrinkageFilter <- function(T, prior.fun=cov.prior.cc, ...)
{
  list(prior.fun=prior.fun, ...)
}

################################### DENOISE ##################################
# p <- create(TawnyPortfolio, h, 90)
# denoise(p,create(SampleFilter))
denoise %when% (estimator %isa% SampleFilter)
denoise %also% (p %isa% TawnyPortfolio)
denoise %as% function(p, estimator)
{
  cov2cor(cov.sample(p$returns))
}

# p <- create(TawnyPortfolio, h, 90)
# denoise(p,create(EmpiricalFilter))
denoise %when% (estimator %isa% EmpiricalFilter)
denoise %also% (p %isa% TawnyPortfolio)
denoise %as% function(p, estimator)
{
  cor.empirical(p$returns)
}

# Customizations:
#   Starting correlation matrix
#   Fit to distribution and get lambda.plus
#   Clean function
#   
# p <- create(TawnyPortfolio, h, 90)
# denoise(p,create(RandomMatrixFilter))
#
# s <- c('FCX','AAPL','JPM','AMZN','VMW','TLT','GLD','FXI','ILF','XOM')
# p <- create(TawnyPortfolio, s)
# w <- rollapply(p, function(x) denoise(x, create(RandomMatrixFilter)))
denoise %when% (estimator %isa% RandomMatrixFilter)
denoise %also% (p %isa% TawnyPortfolio)
denoise %as% function(p, estimator)
{
  #filter.RMT(p$returns, hint=estimator$hint)
  # Use either a fit based on the theoretical shape or use asymptotics for an
  # analytical solution
  #logger(INFO, "Estimating correlation matrix")
  p$correlation <- estimator$cor.fn(p$returns)
  #logger(INFO, "Getting eigenvalues")
  es <- eigen(p$correlation, symmetric=TRUE, only.values=FALSE)
  #logger(INFO, "Applying cutoff")
  lambda.plus <- estimator$cutoff.fn(p$correlation, es, estimator)
  #logger(INFO, "Cleaning matrix")
  estimator$clean.fn(es, lambda.plus)
}

# p <- create(TawnyPortfolio, h, 90)
# denoise(p,create(ShrinkageFilter))
denoise %when% (estimator %isa% ShrinkageFilter)
denoise %also% (p %isa% TawnyPortfolio)
denoise %as% function(p, estimator)
{
  cov2cor(cov_shrink(p$returns, prior.fun=estimator$prior.fun))
}

denoise %when% (estimator %isa% ShrinkageFilter)
denoise %also% (estimator %hasa% market)
denoise %also% (p %isa% TawnyPortfolio)
denoise %as% function(p, estimator)
{
  h <- p$returns
  # TODO: Validation should be in type constructor
  if (!is.null(rownames(h)))
    dates <- as.Date(rownames(h), format='%Y-%m-%d')
  else if (!is.null(index(h)))
    dates <- index(h)
  else
    dates <- names(h)

  # This is not as efficient if the market has yet to be downloaded, but
  # the interface is cleaner. Think about how to balance this better later.
  # TODO: Put in separate function
  if ('character' %in% class(market))
  {
    if (exists(market)) m <- get(market)
    #else m <- getPortfolioReturns(market, obs=obs, start=start, end=end)
    else m <- getPortfolioReturns(market, obs=anylength(h), end=last(dates))
  }
  else m <- market

  if (start(m) > first(dates)) stop("Market data does not span start of h")
  if (end(m) < last(dates)) stop("Market data does not span end of h")

  m <- m[index(m) >= first(dates) & index(m) <= last(dates)]
  if (anylength(m) != anylength(h)) stop("Inconsistent data lengths")

  # Calculate beta
  # TODO: Make this an external portfolio function
  beta.fn <- function(h, m) cov(h, m) / var(m)
  betas <- apply(h, 2, beta.fn, m)

  # Subtract beta * market from returns
  mkt <- matrix(rep(m, ncol(h)), byrow=TRUE, nrow=ncol(h))
  h <- h - t(betas * mkt)


  cov2cor(cov.shrink(h, prior.fun=prior.fun, ...))
}


##------------------------ CORRELATION MATRIX FUNCTIONS ---------------------##
# Clean a correlation matrix based on calculated value of lambda.plus and the
# computed eigenvalues.
# This takes flattened eigenvalues and returns a new cleaned correlation matrix
# Params:
#  es: eigenvalues and vectors
#  lambda.plus: cutoff
#  h: non-normalized returns matrix (only used for labels)
cor.clean <- function(es, lambda.plus=1.6, h=NULL)
{
  e.values <- es$values
  avg <- mean(e.values[e.values < lambda.plus])
  e.values[e.values < lambda.plus] <- avg

  e.vectors <- es$vectors
  c.clean <- e.vectors %*% diag(e.values) %*% t(e.vectors)
  diags <- diag(c.clean) %o% rep(1, nrow(c.clean))
  c.clean <- c.clean / sqrt(diags * t(diags))

  if (! is.null(h) & 'returns' %in% class(h))
  {
    rownames(c.clean) <- anynames(h)
    colnames(c.clean) <- anynames(h)
  }
  c.clean
}

# h TxM zoo returns matrix
# This doesn't subtract the mean (based on the literature)
cor.empirical <- function(h)
{
  # Normalize returns
  ns <- normalize(h)

  # Calculate the correlation matrix
  # E = H H'
  t <- nrow(ns)
  e <- t(ns) %*% ns / t
  class(e) <- c(class(e), 'correlation')
  e
}


# Normalizes a returns matrix such that Var[xit] = 1.
# Assumes TxM, m population, t observations
normalize %when% (h %isa% zoo)
normalize %as% function(h) { apply(h, 2, function(x) x / sd(x)) }
