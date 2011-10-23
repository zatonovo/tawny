create.RandomMatrixEst <- function(T, hint=c(4,1), ...)
{
  list(hint=hint, ...)
}

# Specify market as a symbol if you want to shrink on residuals only
create.ShrinkageEst <- function(T, prior.fun=cov.prior.cc, ...)
{
  list(prior.fun=prior.fun, ...)
}


# p <- create(TawnyPortfolio, h, 90)
# denoise(p,create(SampleEst))
denoise.sample %when% (estimator %isa% SampleEst & p %isa% TawnyPortfolio)
denoise.sample <- function(p, estimator)
{
  cov2cor(cov.sample(p$returns))
}

# p <- create(TawnyPortfolio, h, 90)
# denoise(p,create(EmpiricalEst))
denoise.emp %when% (estimator %isa% EmpiricalEst & p %isa% TawnyPortfolio)
denoise.emp <- function(p, estimator)
{
  cor.empirical(p$returns)
}

# p <- create(TawnyPortfolio, h, 90)
# denoise(p,create(RandomMatrixEst))
denoise.rmt %when% (estimator %isa% RandomMatrixEst & p %isa% TawnyPortfolio)
denoise.rmt <- function(p, estimator)
{
  filter.RMT(p$returns, hint=estimator$hint)
  # Use either a fit based on the theoretical shape or use asymptotics for an
  # analytical solution
}

# p <- create(TawnyPortfolio, h, 90)
# denoise(p,create(ShrinkageEst))
denoise.shrink %when% (estimator %isa% ShrinkageEst & p %isa% TawnyPortfolio)
denoise.shrink <- function(p, estimator)
{
  cov2cor(cov.shrink(p$returns, prior.fun=estimator$prior.fun))
}

denoise.shrinkm %when% (estimator %isa% ShrinkageEst & estimator %hasa% market)
denoise.shrinkm <- function(p, estimator)
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

