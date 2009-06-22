library(zoo)
library(quantmod)

# No filtering
getCorFilter.Sample <- function() { function(h) cov.sample(h) }

# This acts as a control case with no cleaning
getCorFilter.Raw <- function()
{
  function(h) return(cor.empirical(h))
}

# This is the default filter
getCorFilter.RMT <- function(hint=c(4,1), ...)
{
  # h is a zoo object and will be t x m
  #function(h) { filter.RMT(t(h), breaks, hint) }
  function(h) return(filter.RMT(h, hint=hint, ...))
}

# Return a correlation matrix generator that is compatible with the portfolio
# optimizer
# Example
#   ws <- optimizePortfolio(ys, 100, getCorFilter.Shrinkage())
#   plotPerformance(ys,ws)
getCorFilter.Shrinkage <- function(prior.fun=cov.prior.cc, ...)
{
  function(h) return(cov2cor(cov.shrink(h, prior.fun=prior.fun, ...)))
}

# Use shrinkage but subtract the beta * market for each stock and look at
# residuals only
getCorFilter.ShrinkageM <- function(market, prior.fun=cov.prior.cc, ...)
{

  fn <- function(h)
  {
    dates <- as.Date(rownames(h), format='%Y-%m-%d')
    # This is not as efficient if the market has yet to be downloaded, but
    # the interface is cleaner. Think about how to balance this better later.
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
    beta.fn <- function(h, m) cov(h, m) / var(m)
    betas <- apply(h, 2, beta.fn, m)

    # Subtract beta * market from returns
    mkt <- matrix(rep(m, ncol(h)), byrow=TRUE, nrow=ncol(h))
    h <- h - t(betas * mkt)


    cov2cor(cov.shrink(h, prior.fun=prior.fun, ...))
  }
  fn
}

