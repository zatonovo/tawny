# High level functions
# tny.optimize(Portfolio, Model)
# denoise(Portfolio, Model)
# backtest(Portfolio, Denoiser, Optimizer)


# h <- getPortfolioReturns(getIndexComposition('^DJI'), obs=400, reload=TRUE)
# ws <- optimizePortfolio(h, 350, RandomMatrixDenoiser() )
# pf <- plotPerformance(h, ws, 350, y.min=-0.4)
# ef <- compare.EqualWeighted(h, 350, y.min=-0.4)
# mf <- compare.Market('^GSPC',400,350, y.min=-0.4)



# Optimizes a returns series over a window.
# s <- c('FCX','AAPL','JPM','AMZN','VMW','TLT','GLD','FXI','ILF','XOM')
# p <- TawnyPortfolio(s)
# ws <- optimizePortfolio(p, RandomMatrixDenoiser())

optimizePortfolio(p, estimator, optimizer) %::% TawnyPortfolio : a : Function : z
optimizePortfolio(p, estimator, optimizer=p.optimize) %as% {
  my.optimizer <- function(p)
  {
    flog.debug("Getting correlation matrix")
    #cor.mat <- cor.gen(h.window, ...)
    cor.mat <- denoise(p, estimator)

    flog.debug("Optimizing portfolio")
    optimizer(p$returns, cor.mat)
  }

  flog.info("Optimizing portfolio for [%s, %s]",
         format(start(p)), format(end(p)) )
  ws <- rollapply(p, my.optimizer)
  xts(ws, order.by=as.Date(rownames(ws)))
}

# This is for backwards compatibility
optimizePortfolio(h, window, estimator, ...) %as%
{
  flog.info("Note: This interface is deprecated")
  p <- TawnyPortfolio(h, window)
  optimizePortfolio(p, estimator)
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
    flog.warn("Unable to invert correlation matrix. Returning zeros.")
    return(rep(0, ncol(h)))
  }

  sdna <- function(x) { sd(x,na.rm=TRUE) }
  h.vol <- apply(h, 2, sdna)
  sig.clean <- c.inv / (h.vol[col(c.inv)] * h.vol[row(c.inv)])
  norm <- sum(sig.clean)

  apply(sig.clean,1,sum) / norm
}


