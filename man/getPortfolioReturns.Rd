\name{getPortfolioReturns}
\alias{getIndexComposition}
\alias{getPortfolioReturns}
\alias{ensure}
\title{ Utility functions for creating portfolios of returns and other functions }
\description{
  Gets portfolio returns from closing prices (configurable). This uses quantmod
  under the hood to retrieve prices and construct returns based on a 
  configurable transform.

  Also included is a function that returns the composition of select indexes
  that can be used in conjunction with getPortfolioReturns() to get the underlying
  returns of the given index.

  Additionally, there is a utility function that ensures that symbols have been
  properly loaded.
}
\usage{
getIndexComposition(ticker = '^GSPC', hint = NA, src = 'yahoo')
getPortfolioReturns(symbols, obs = NULL, start = NULL, end = Sys.Date(), fun = function(x) Delt(Cl(x)), reload = FALSE, na.value = NA, ...)
ensure(serie, src = 'FRED', reload = FALSE, ...)
}
\arguments{
  \item{ticker}{ The ticker of the index. For best mileage, use Yahoo! compatible tickers (including the caret prefix). }
  \item{hint}{ A hint that specifies the number of assets in the index. If omitted, a default will be used based on pre-configured data for common indexes. }
  \item{src}{ The data vendor to use. Defaults to yahoo but could work with google }
  \item{symbols}{ A vector (or scalar) of tickers to download }
  \item{obs}{ The number of observations to get }
  \item{start}{ Alternatively, a start date can be used to specify the beginning of a range to download. }
  \item{end}{ The end date of the range. Defaults to current date. }
  \item{fun}{ A transform applied to the downloaded data. The default is to calculate returns on the close. }
  \item{reload}{ Whether to reload data or just download missing data }
  \item{na.value}{ What value to use if the resulting portfolio has NAs. The default is to omit any assets containing NA values. }
  \item{serie}{ A vector (or scalar) of tickers to ensure exist in the current environment }
  \item{\dots}{ Additional parameters to pass to getSymbols }
}
\details{
  Typically only getPortfolioReturns and getIndexComposition will be used on a 
  regular basis.
  
  Ensure isn't as useful as initially conceived given that naming collisions
  have caused numerous issues. The code now uses auto.assign=FALSE in the 
  underlying getSymbols call.
}
\value{
  getPortfolioReturns returns a TxM xts object of asset returns.

  getIndexComposition returns a vector of asset symbols (i.e. tickers).

  ensure returns nothing.
}
\author{ Brian Lee Yung Rowe}
\seealso{ \code{\link{tawny}} }
\examples{
# This is autorun outside of examples
tawny:::.init()

# Get a portfolio
h <- getPortfolioReturns(c('A','AA','AAPL'), obs=150)

# Get an index portfolio
h <- getPortfolioReturns(getIndexComposition('^DJI'), obs=100, reload=TRUE)

# Doesn't work because of numerical symbols - need to fix
#h <- getPortfolioReturns(getIndexComposition('^HSI'), obs=100, reload=TRUE)

# Ensure that some assets exist
ensure(c('K','JNPR'), src='yahoo')
}
\keyword{ ts }
