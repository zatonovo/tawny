# High level functions
# tny.optimize(Portfolio, Model)
# denoise(Portfolio, Model)
# backtest(Portfolio, Denoiser, Optimizer)


# h <- getPortfolioReturns(getIndexComposition('^DJI'), obs=400, reload=TRUE)
# ws <- optimizePortfolio(h, 350, getCorFilter.RMT() )
# pf <- plotPerformance(h, ws, 350, y.min=-0.4)
# ef <- compare.EqualWeighted(h, 350, y.min=-0.4)
# mf <- compare.Market('^GSPC',400,350, y.min=-0.4)


# Example
#  h <- getPortfolioReturns(c('MER','C','GS','MS','BAC','AAPL','GOOG','MSFT','IBM','CSCO'),150)
## Not run:
#  h <- getPortfolioReturns(c('MER','C','GS','MS','BAC','WFC',
#    'ORCL','YHOO','T','AAPL','GOOG','MSFT','IBM','CSCO','INTC',
#    'GM','F','CAT','MMM','XOM','CVX','RIG','USO','ABX',
#    'AMGN','HGSI','PFE','JNJ', 'FSLR','STP'), 150, reload=TRUE)
#  c.gen <- getCorFilter.RMT(hint=c(4,1))
#  weights <- optimizePortfolio(h, 100, c.gen)
## End(Not run)
# NOTE: For zoo compatibility, need to use a t x m matrix for h
# Attempls to dispatch based on characteristics of h.

# Optimizes a returns series over a window.
# s <- c('FCX','AAPL','JPM','AMZN','VMW','TLT','GLD','FXI','ILF','XOM')
# p <- create(TawnyPortfolio, s)
# ws <- optimizePortfolio(p, create(RandomMatrixFilter))
optimizePortfolio %when% (p %isa% TawnyPortfolio)
optimizePortfolio %as% function(p, estimator)
{
  optimizePortfolio(p, estimator, p.optimize)
}

optimizePortfolio %when% (p %isa% TawnyPortfolio)
optimizePortfolio %as% function(p, estimator, optimizer)
{
  my.optimizer <- function(p)
  {
    logger(DEBUG, "Getting correlation matrix")
    #cor.mat <- cor.gen(h.window, ...)
    cor.mat <- denoise(p, estimator)

    logger(DEBUG, "Optimizing portfolio")
    optimizer(p$returns, cor.mat)
  }

  logger(INFO, sprintf("Optimizing portfolio for %s-%s",
         format(start(p)), format(end(p)) ))
  ws <- rollapply(p, my.optimizer)
  xts(ws, order.by=as.Date(rownames(ws)))
}

# This is for backwards compatibility
#optimizePortfolio <- function(h, window, cor.gen=getCorFilter.RMT(), ...)
optimizePortfolio %when% TRUE
optimizePortfolio %as% function(h, window, estimator, ...)
{
  if (! 'zoo' %in% class(h))
  { cat("WARNING: Zoo objects are preferred for dimensional safety\n") }

  my.optimizer <- function(h.window)
  {
    logger(DEBUG, "Getting correlation matrix")
    #cor.mat <- cor.gen(h.window, ...)
    cor.mat <- denoise(h, estimator, ...)

    logger(DEBUG, "Optimizing portfolio")
    p.optimize(h.window, cor.mat)
  }

  logger(INFO, sprintf("Optimizing portfolio for %s-%s",format(start(h)),format(end(h))))
  ws <- rollapply(h, window, my.optimizer, by.column=FALSE, align='right')
  xts(ws, order.by=index(ws))
}




##------------------------ PRIVATE PORTFOLIO FUNCTIONS ----------------------##
# Optimize the portfolio by identifying the minimum variance portfolio
# Params:
#  h: TxM zoo returns matrix
#  c.denoised: Denoised correlation matrix (m x m)
# Returns:
#  Optimized weights vector for portfolio holdings. Sum of weights = 1.
p.optimize <- function(h, c.denoised)
{
  #if (logLevel() > 0) cat("Optimizing for",as.Date(end(h)),"\n")

  try(c.inv <- solve(c.denoised))
  if (!exists('c.inv'))
  {
    logger(WARN, "Unable to invert correlation matrix. Returning zeros.")
    return(rep(0, ncol(h)))
  }

  sdna <- function(x) { sd(x,na.rm=TRUE) }
  h.vol <- apply(h, 2, sdna)
  sig.clean <- c.inv / (h.vol[col(c.inv)] * h.vol[row(c.inv)])
  norm <- sum(sig.clean)

  apply(sig.clean,1,sum) / norm
}


