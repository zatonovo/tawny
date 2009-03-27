\name{optimizePortfolio}
\alias{optimizePortfolio}
\alias{optimizePortfolio.RMT}
\alias{getCorFilter.RMT}
\alias{p.optimize}
\title{ Optimize a portfolio using the specified correlation matrix generator }
\description{
  Performs basic minimum variance optimization on the given portfolio over the 
  given window returning a zoo object of portfolio weights. The window must be 
  less than or equal to the number of observations in h (s.t. if they are equal
  then only one resultant weight vector will be returned).
}
\usage{
optimizePortfolio(h, window, cor.gen, ...)
optimizePortfolio.RMT(h, window = NULL, breaks = NULL, hint = c(4, 1))
}
\arguments{
  \item{h}{ A zoo object representing a portfolio with dimensions T x M }
  \item{window}{ The size of the window to use when generating a correlation matrix }
  \item{cor.gen}{ A function that takes a portfolio with window observations and generates a correlation matrix }
  \item{breaks}{ The explicit breaks to use when generating the histogram }
  \item{hint}{ A hint for the correct value of Q and sigma }
  \item{\dots}{ Additional parameters to pass to the cor.gen function }
}
\details{
  This is the primary entry point to using the tawny package. This function 
  calculates the portfolio weights over the portfolio based on a rolling window.
  Given M assets in the portfolio, T total observations, and a window of length
  t, the resulting weights object will have dimensions T - t + 1 x M.

  The weights matrix can then be analyzed to calculate standard portfolio
  performance metrics. A simple analytics function is provided so that 
  cumulative returns can be easily viewed, although for more sophisticated
  analysis other packages should be used.

  In theory any compatible correlation matrix generator can be used (and has in
  practice to test against proprietary risk models) and the function will
  generate portfolio weights accordingly. To leverage the remainder of the
  package, the \code{getCorFilter.RMT} function or \code{getCorFilter.Shrinkage}
  should be called. These wrappers are somewhat superfluous but does provide 
  some utility by ensuring compatibility with the underlying RMT code that uses
  transposed matrices (pre-zoo integration). Additionally, by way of closures
  these functions are used to store hints to the optimizer and any final data
  massaging, potentially cleaning up code but admittedly can be serviced via 
  the normal dots mechanism.
  In the future, the default will be a direct handle to the underlying function
  once the rest of the code is converted to zoo.

  The secondary function \code{optimizePortfolio.RMT} exists to optimize the 
  correlation matrix using RMT exclusively. This is a more direct route to 
  accessing the RMT functionality and might be more convenient to use. The 
  intention is that the base optimizePortfolio function becomes a generic 
  function that passes on to specific implementations, but the mechanics haven't
  been worked out yet. It is also possible to extract the optimizer and pass 
  that in explicitly as a function.
}
\value{
  A weights zoo object with T - t + 1 dates and M assets. The dates are aligned
  to the end date.
}
\author{ Brian Lee Yung Rowe}
\examples{
data(sp500.subset)
ys <- sp500.subset
ws <- optimizePortfolio(ys, 150, getCorFilter.RMT() )
}
\keyword{ ts }
