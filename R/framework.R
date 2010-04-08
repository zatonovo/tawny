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
optimizePortfolio <- function(h, window, cor.gen=getCorFilter.RMT(), ...)
{
  if (! 'zoo' %in% class(h))
  { cat("WARNING: Zoo objects are preferred for dimensional safety\n") }

  my.optimizer <- function(h.window)
  {
    logger(DEBUG, "Getting correlation matrix")
    cor.mat <- cor.gen(h.window, ...)

    logger(DEBUG, "Optimizing portfolio")
    p.optimize(h.window, cor.mat)
  }

  logger(INFO, sprintf("Optimizing portfolio for %s-%s",format(start(h)),format(end(h))))
  ws <- rollapply(h, window, my.optimizer, by.column=FALSE, align='right')
  xts(ws, order.by=index(ws))
}

##----------------------------- TESTING FUNCTIONS ---------------------------##
# Plot the performance of the portfolio
# Params:
#  h: A returns matrix. This can be either in-sample or out-of-sample
#  weights: Either a vector or an m x t matrix of weights
# Usage:
#  h <- some.returns.matrix # m x t
#  ws <- optimizePortfolio.RMT(h)
#  plotPerformance(h,ws)
# Example:
#  load('RandomMatrix/spxReturnsQ5.rData')
#  spx.small <- t(as.matrix(spx.retq5[1:250,][,1:50]))
#  spx.small.history <- t(spx.retq5)[1:50,]
# NOTE: for zoo compatibility,
#  h is t x m
#  w is t x m
# ws <- optimizePortfolio(spx[1:350, 1:50], 250, getCorFilter.RMT())
# ps <- plotPerformance(spx[1:350, 1:50], ws)

# perf <- plotPerformance(h,weights,100, y.min=-0.75)
# The weights are currently expected to be the same as the last date in the
# window. Consequently, the weights necessarily are valid for T+1 at the 
# earliest.
plotPerformance <- function(h, weights, window=NULL, rf.rate=0.01,
  new.plot=TRUE, y.min=-0.25, y.max=0.25, bg=NULL,
  name='', color='red', colors=c(), legend.loc='bottomright',...)
{
  if (is.null(window)) { window <- anylength(h) - anylength(weights) + 1 }

  if (! 'zoo' %in% class(h))
  { logger(WARN,"Zoo objects are preferred for dimensional safety") }

  if (anylength(h) != anylength(weights) + window - 1)
  { logger(WARN,"Dimensions are inconsistent. This can cause errors") }

  ts.rets <- portfolioReturns(h, weights)
  stats <- portfolioPerformance(ts.rets, rf.rate)

  # Plot this output
  logger(TRACE, sprintf("y range = [%s,%s]",y.min,y.max))
  #cat("x range = [",as.Date(start(xaxis)),",",as.Date(end(xaxis)),"]\n")

  #yrange=c(min(y.min, ts.perf-0.05), max(y.max, ts.perf+0.05))
  # Need to synchronize y bounds for multiple charts
  yrange=c(y.min, y.max)
  #if (! is.null(bg)) par(bg=bg, lab=c(7,7,7))
  #else par(lab=c(7,7,7))

  logger(TRACE, "Plotting chart")
  if (anylength(dev.list()) > 0 & new.plot & is.null(bg))
  {
    par(lab=c(7,7,7), new=new.plot)
    lines(stats$cum.returns, col=color, ...)
  }
  else
  {
    if (!is.null(bg)) par(bg=bg)
    par(lab=c(7,7,7))
    plot(stats$cum.returns,
      type="l",col=color,ylim=yrange,ylab="Return",xlab="Time", ...);
  }

  # Ignore line types for now
  #stats$ltys <- c(ltys, lty)
  stats$colors <- c(colors, name=color)
  names(stats$colors)[anylength(stats$colors)] <- name
  grid()
  #legend('topright', legend=names(stats$colors), lwd=2, cex=1.0, 
  #  col=stats$colors, border.col='slategray4', inset=0.02)
  legend(legend.loc, legend=names(stats$colors), lwd=2, cex=1.0, 
    col=stats$colors, inset=0.02)

  stats
}

# Calculate portfolio returns based
portfolioReturns <- function(h, weights)
{
  logger(TRACE, sprintf("class(index(h)):%s",class(index(h))))

  # Shift dates so weights are used on following date's data for out-of-sample
  # performance
  w.index <- c(index(weights[2:anylength(weights)]), end(weights) + 1)
  index(weights) <- w.index

  h.trim <- h[index(h) %in% index(weights)]
  ts.rets <- apply(zoo(h.trim) * zoo(weights), 1, sum)

  # This is in here to fix some strange behavior related to rownames vs index
  # in zoo objects and how they are used after an apply function
  #names(ts.rets) <- index(h.trim)
  #ts.rets <- zoo(ts.rets, order.by=as.Date(names(ts.rets)))
  ts.rets <- zoo(ts.rets, order.by=index(h.trim))

  # This causes problems
  #ts.rets <- xts(t(h.trim * weights) %*% rep(1, anylength(h.trim)), order.by=index(h.trim))
  if (any(is.na(ts.rets)))
  {
    logger(WARN,"Filling NA returns with 0")
    ts.rets[is.na(ts.rets)] <- 0
  }

  logger(DEBUG, sprintf("Returns count:%s",anylength(ts.rets)))

  return(ts.rets)
}

# p - portfolio returns
portfolioPerformance <- function(p, rf.rate=0.01)
{
  #require(PerformanceAnalytics, quietly=TRUE)
  #ts.perf <- cumprod(1+p) - 1
  #xaxis <- as.Date(index(ts.perf))
  #xaxis <- as.Date(index(p))
  days <- length(p)

  stats <- list()
  stats$daily.returns <- p
  stats$cum.returns <- cumprod(1 + stats$daily.returns) - 1
  stats$period.stdev <- as.numeric(sd(p, na.rm=TRUE))
  stats$annual.stdev <- as.numeric(sd(p, na.rm=TRUE) * 252^0.5)
  #stats$drawdown <- maxDrawdown(p)
  stats$avg.return <- mean(stats$daily.returns)
  stats$period.return <- as.numeric(last(stats$cum.returns))
  stats$annual.return <- (1 + stats$period.return)^(252/days) - 1
  stats$sharpe.ratio <- (stats$annual.return - rf.rate) / stats$annual.stdev
  #prf <- (1 + rf.rate)^(window/252) - 1
  #stats$sharpe.ratio <- SharpeRatio.annualized(stats$daily.returns, rf=prf, scale=252)

  return(stats)
}

# Compare with the given market
# mkt.perf <- compare.Market('^DJI',150,100)
# mkt.perf <- compare.Market(DJI,150,100, y.min=-0.75)
compare.Market <- function(market, obs, window, end=Sys.Date(), color='#44bc43', ...)
{
  if (is.character(market))
  {
    start <- end - (10 + obs * 365/250)
    mkt <- getSymbols(market, src='yahoo',from=start,to=end, auto.assign=FALSE)
  }
  else mkt <- market

  mkt.ret <- Delt(Cl(mkt))
  mkt.ret <- mkt.ret[index(mkt.ret) <= end]
  mkt.ret <- tail(mkt.ret, obs)

  w.count <- obs - window + 1
  weights <-
    xts(matrix(1,ncol=1,nrow=w.count), order.by=index(tail(mkt.ret, w.count)))

  plotPerformance(mkt.ret, weights, window, color=color, name='market', ...)
}

# Compare the performance with an equal weighted portfolio
# eq.perf <- compare.EqualWeighted(h, 100, y.min=-0.75)
compare.EqualWeighted <- function(h, window, color='#342a31', ...)
{
  my.weights <- function(x) rep(1/ncol(x), ncol(x))
  weights <- xts(rollapply(h, window, my.weights, by.column=FALSE, align='right'), order.by=index(h[window:anylength(h)]))
  logger(INFO, sprintf("weights: [%s,%s]",as.Date(start(weights)),as.Date(end(weights))))

  #index(weights) <- index(h[window:anylength(h)])
  plotPerformance(h, weights, window, color=color, name='naive', ...)
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


