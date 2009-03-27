# Look at rmetrics.org
#library(PerformanceAnalytics)
library(futile)
library(zoo)
library(quantmod)

# TODO
# x. Modify optimizePortfolio to use rollapply
# x. Modify optimizePortfolio.RMT to use rollapply
# x. Create portfolio function using zoo/xts
# x. Auto generate equal-weighted portfolio to compare with
# . Time-awareness in ensure so that keeping a session open does the right
#   thing
# . Verify/stress test code
# x. Clean up code and comments
# x. Clean up logging
# . Provide example using quantmod
# x. Use PerformanceAnalytics for portfolio performance (LATER)
# x. Possibly modify chart (filter.RMT, mp.theory) to use quantmod -
#   doesn't seem to be compatible given non-time series nature of some of the 
#   plots
# x. Disable charting of fit as an option
# x. Dates for weights are incorrect (why?)
# x. Line up dates

# h <- getPortfolioReturns(getIndexComposition('^DJI'), obs=400, reload=TRUE)
# ws <- optimizePortfolio(h, 350, getCorFilter.RMT() )
# pf <- plotPerformance(h, ws, 350, y.min=-0.4)
# ef <- compare.EqualWeighted(h, 350, y.min=-0.4)
# mf <- compare.Market('^GSPC',400,350, y.min=-0.4)



##---------------------------- PUBLIC FUNCTIONS -----------------------------##
# Optimize a portfolio using the specified portfolio generator and correlation
# matrix generator. This is a generic optimizer that allows any custom generator
# of correlation matrices to be used.

# Returns current log level of package
#logLevel <- function(new.level=NULL)
#{
#  if (! is.null(new.level)) { options('log.level'=new.level) }
#
#  if (is.null(getOption('log.level'))) { return(0) }
#  return(getOption('log.level'))
#}

# Set with options('use.plots'=FALSE)
# Defaults to TRUE
#usePlots <- function(new.val=NULL)
#{
#  if (! is.null(new.val)) { options('use.plots'=new.val) }
#
#  if (is.null(getOption('use.plots'))) { return(TRUE) }
#  return(getOption('use.plots'))
#}

# Ensures that a given series exists and downloads via quantmod if it doesn't
# serie - Either a string or a list
ensure <- function(serie, src='FRED', reload=FALSE, ...)
{
  for (series in serie)
  {
    if (! reload)
    {
      cleaned <- sub('^','', series, fixed=TRUE)
      if (exists(cleaned)) { next }
      cat("Loading new symbol",series,"from",src,"\n")
    }
    else
    {
      cat("(Re)loading symbol",series,"from",src,"\n")
    }
    try(getSymbols(series, src=src, ...))
  }
}

# Example
# Get SP500 components
#   sp500.idx <- getIndexComposition()
# Get DOW components
#   dow.idx <- getIndexComposition('^DJI')
# Get FTSE components
#   ftse.idx <- getIndexComposition('^FTSE')
# Get HSI components
#   hsi.idx <- getIndexComposition('^HSI')
# h <- getPortfolioReturns(getIndexComposition('^DJI'), obs=100)
getIndexComposition <- function(ticker='^GSPC', hint=NA, src='yahoo')
{
  if (is.na(hint))
  {
    hints <- c(500, 30, 102, 42)
    names(hints) <- c('^GSPC', '^DJI', '^FTSE', '^HSI')
    hint <- hints[ticker]
  }

  # http://download.finance.yahoo.com/d/quotes.csv?s=@%5EGSPC&f=sl1d1t1c1ohgv&e=.csv&h=0
  base <- 'http://download.finance.yahoo.com/d/quotes.csv?s=@'
  formats <- '&f=sl1d1t1c1ohgv&e=.csv&h='

  comp <- NULL
  pages = max(1, hint %/% 50)
  for (page in 1:pages)
  {
    start <- (page-1) * 50 + 1
    url <- paste(base, ticker, formats, start, sep='')
    cat("Loading page",page,"for",ticker,"\n")
    data <- read.csv(url, header=FALSE)

    # This is here due to a bug in Yahoo's download where the first record gets
    # duplicated in each subsequent page
    idx = 2; if (page == 1) { idx = 1 }
    comp <- rbind(comp, data[idx:anylength(data),])

  }
  as.character(comp[,1])
}

# This produces a portfolio in matrix format (t x m) as a zoo class. 
# Params
#  symbols: A vector of symbols to retrieve. This uses quantmod to retrieve
#    the data.
#  obs: The number of observations that you want. Use this if you want the 
#    number of points to be explicit. Either obs or start is required.
#  start: The start date, if you know that explicitly. Using this will ensure
#    that the data points are bound to the given range but the precise number
#    of points will be determined by the number of trading days.
#  end: The most recent date of observation. Defaults to current day.
#  fun: A function to use on each symbol time series. Defaults to Cl to operate
#    on close data. For expected behavior, your function should only return
#    one time series.
# TODO: 
#  Fix names
#  Add method to add other portfolio elements (such as synthetic securities)
getPortfolioReturns <- function(symbols, obs=NULL, start=NULL, end=Sys.Date(),
  fun=function(x) Delt(Cl(x)), reload=FALSE, na.value=NA, ...)
{
  if (is.null(start) & is.null(obs)) { stop("Either obs or start must be set") }
  end <- as.Date(end)

  # Estimate calendar days from windowed business days. The 10 is there to
  # ensure enough points, which get trimmed later
  if (is.null(start)) { start <- end - (10 + obs * 365/250) }

  # Load symbols - The problem with this is that it is not time-aware. Need to
  # think about what to do with it.
  ensure(symbols, src='yahoo', reload=reload, from=start, to=end, ...)

  # Merge into a single zoo object
  p <- xts(order.by=last(index(get(symbols[1]))))
  for (s in symbols)
  {
    if (!exists(s)) next
    if (logLevel() > 0) cat("Binding",s,".")
    raw <- fun(get(s))
    if (logLevel() > 0) cat("Start:",start(raw),"; end:",end(raw),"\n")
    a <- xts(raw, order.by=index(get(s)))
    p <- cbind(p, a[2:anylength(a)])
  }
  anynames(p) <- symbols
  # First remove dates that have primarily NAs (probably bad data)
  o.dates <- rownames(p)
  p <- p[apply(p, 1, function(x) sum(x, na.rm=TRUE) != 0), ]
  cat("Removed suspected bad dates",setdiff(o.dates,rownames(p)),"\n")

  if (! is.na(na.value))
  {
    #for (s in symbols) p[,s][is.na(p[,s])] <- na.value
    p[is.na(p)] <- 0
    cat("Replaced NAs with",na.value,"\n")
  }
  else
  {
    # NOTE: This has consistency issues when comparing with a market index
    o.dates <- rownames(p)
    p <- p[apply(p, 1, function(x) sum(is.na(x)) < 0.1 * length(x) ), ]
    cat("Removed dates with too many NAs",setdiff(o.dates,rownames(p)),"\n")

    # Now remove columns with NAs
    nas <- apply(p, 2, function(x) !any(is.na(x)) )
    p <- p[,which(nas == TRUE)]
    cat("Removed symbols with NAs:",setdiff(symbols,anynames(p)),"\n")
  }

  if (is.null(obs)) { return(p[paste(start,end, sep='::')]) }

  p <- p[index(p) <= end]
  idx.inf <- anylength(p) - min(anylength(p), obs) + 1
  idx.sup <- anylength(p)
  if (logLevel() > 0) cat("Returning rows [",idx.inf,",",idx.sup,"]\n")
  
  cat("Loaded portfolio with",ncol(p),"assets\n")
  p[idx.inf:idx.sup, ]
}

# Return a portion of a matrix. This is useful for debugging.
peek <- function(x, upper=5, lower=1)
{
  return(x[lower:upper,lower:upper])
}

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
optimizePortfolio <- function(h, window, cor.gen=getCorFilter.RMT(), ...)
{
  if (! 'zoo' %in% class(h))
  { cat("WARNING: Zoo objects are preferred for dimensional safety\n") }

  log.level <- logLevel()
  my.optimizer <- function(h.window)
  {
    if (log.level > 2) { cat("Getting correlation matrix\n") }
    cor.mat <- cor.gen(h.window, ...)

    if (log.level > 2) { cat("Optimizing portfolio\n") }
    p.optimize(h.window, cor.mat)
  }

  if (log.level > 0)
  {
    cat("Optimizing portfolio for",format(start(h)),"-",format(end(h)),"\n")
  }
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
  new.plot=TRUE, y.min=-.25, y.max=.25, bg=NULL,
  name='', color='red', colors=c(), lty=1, ltys=c(), ...)
{
  require(PerformanceAnalytics, quietly=TRUE)
  if (is.null(window)) { window <- anylength(h) - anylength(weights) }

  if (! 'zoo' %in% class(h))
  { cat("WARNING: Zoo objects are preferred for dimensional safety\n") }

  if (anylength(h) != anylength(weights) + window - 1)
  { cat("WARNING: Dimensions are inconsistent. This can cause errors\n") }

  log.level <- logLevel()

  if (log.level > 5) { cat("class(index(h)):",class(index(h)),"\n") }

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
    cat("WARNING: Filling NA returns with 0\n")
    ts.rets[is.na(ts.rets)] <- 0
  }

  #ts.perf <- cumprod(1+ts.rets) - 1
  #xaxis <- as.Date(index(ts.perf))
  #xaxis <- as.Date(index(ts.rets))

  if (log.level > 2)
  {
    cat("Returns count:",anylength(ts.rets),"\n")
    cat("Date count:",anylength(xaxis),"\n")
  }

  stats <- list()
  stats$daily.returns <- ts.rets
  stats$cum.returns <- cumprod(1 + stats$daily.returns) - 1
  stats$period.stdev <- as.numeric(sd(ts.rets, na.rm=TRUE))
  stats$annual.stdev <- as.numeric(sd(ts.rets, na.rm=TRUE) * (252/window)^0.5)
  #stats$drawdown <- maxDrawdown(ts.rets)
  stats$avg.return <- mean(stats$daily.returns)
  stats$period.return <- as.numeric(last(stats$cum.returns))
  stats$annual.return <- (1 + stats$period.return)^(252/window) - 1
  #stats$sharpe.ratio <- (stats$annual.return - rf.rate) / stats$annual.stdev
  prf <- (1 + rf.rate)^(window/252) - 1
  stats$sharpe.ratio <- SharpeRatio.annualized(stats$daily.returns, rf=prf, scale=252)

  # Plot this output
  if (log.level > 3)
  {
    cat("y range = [",y.min,",",y.max,"]\n")
    cat("x range = [",as.Date(start(xaxis)),",",as.Date(end(xaxis)),"]\n")
  }
  #yrange=c(min(y.min, ts.perf-0.05), max(y.max, ts.perf+0.05))
  # Need to synchronize y bounds for multiple charts
  yrange=c(y.min, y.max)
  if (! is.null(bg)) par(bg=bg, lab=c(7,7,7))
  else par(lab=c(7,7,7))

  if (log.level > 3) { cat("Plotting chart\n") }
  #plot(xaxis, stats$cum.returns,
  if (anylength(dev.list()) > 0 & new.plot)
  {
    par(new=new.plot)
    lines(stats$cum.returns, col=color, ...)
  }
  else
  {
    plot(stats$cum.returns,
      type="l",col=color,ylim=yrange,ylab="Return",xlab="Time", ...);
  }

  stats$ltys <- c(ltys, lty)
  stats$colors <- c(colors, name=color)
  names(stats$colors)[anylength(stats$colors)] <- name
  grid()
  legend('topright', legend=names(stats$colors), lwd=2, cex=1.0, 
    col=stats$colors, border.col='slategray4', inset=0.02)
  #legend('topright', legend=names(stats$colors), lwd=1, cex=1.0, 
  #  col=stats$colors, lty=stats$ltys)

  stats
}

# Compare with the given market
# mkt.perf <- compare.Market('^DJI',150,100)
# mkt.perf <- compare.Market(DJI,150,100, y.min=-0.75)
compare.Market <- function(market, obs, window, end=Sys.Date(), color='#44bc43', ...)
{
  if (is.character(market))
  {
    start <- end - (10 + obs * 365/250)
    ensure(market, src='yahoo', from=start, to=end)
    mkt <- get(sub('^','', market, fixed=TRUE))
  }
  else mkt <- market
  #mkt.trim <- mkt[index(mkt) <= end,]
  #mkt.trim <- mkt.trim[(anylength(mkt.trim)-obs+1) : anylength(mkt.trim)]

  mkt.ret <- Delt(Cl(mkt))
  mkt.ret <- mkt.ret[index(mkt.ret) <= end]
  mkt.ret <- mkt.ret[(anylength(mkt.ret)-(obs-window+1)):anylength(mkt.ret)]

  weights <- xts(rep(1,anylength(mkt.ret)), order.by=index(mkt.ret))
  # This is just a cheap way to get an xts object with all 1s
  #weights <- Hi(mkt.trim[window:anylength(mkt.ret)]) ^ 0
  #index(weights) <- index(mkt.ret[window:anylength(mkt.ret)])
  #weights <- xts(rep(1, anylength(mkt.trim)-obs), order.by=index(mkt.ret[window:anylength(mkt.ret)]))

  plotPerformance(mkt.ret, weights, window, color=color, name='market', ...)
}

# Compare the performance with an equal weighted portfolio
# eq.perf <- compare.EqualWeighted(h, 100, y.min=-0.75)
compare.EqualWeighted <- function(h, window, color='#342a31', ...)
{
  my.weights <- function(x) rep(1/ncol(x), ncol(x))
  weights <- xts(rollapply(h, window, my.weights, by.column=FALSE, align='right'), order.by=index(h[window:anylength(h)]))
  if (logLevel() > 0) 
  { cat("weights: [",as.Date(start(weights)),",",as.Date(end(weights)),"]\n") }

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
    cat("Unable to invert correlation matrix. Returning zeros.\n")
    return(rep(0, ncol(h)))
  }

  sdna <- function(x) { sd(x,na.rm=TRUE) }
  h.vol <- apply(h, 2, sdna)
  sig.clean <- c.inv / (h.vol[col(c.inv)] * h.vol[row(c.inv)])
  norm <- sum(sig.clean)

  apply(sig.clean,1,sum) / norm
}


