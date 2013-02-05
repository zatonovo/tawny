\name{cov_shrink}
\alias{cov_shrink}
\alias{cov_sample}
\alias{cov.prior.cc}
\alias{cov.prior.identity}
\alias{cor.mean}
\alias{shrinkage.p}
\alias{shrinkage.r}
\alias{shrinkage.c}
\alias{shrinkage.intensity}
\alias{cov.shrink} % Deprecated
\alias{cov.sample} % Deprecated
\title{ Shrink the covariance matrix towards some global mean }
\description{
  This performs a covariance shrinkage estimation as specified in Ledoit and
  Wolf. Using within the larger framework only requires using the 
  ShrinkageDenoiser type, which handles the work of constructing a
  shrinkage estimate of the covariance matrix of returns (and consequently its
  corresponding correlation matrix).
}
\usage{
cov_shrink(...)

cov_sample(...)

cov.prior.cc(S)

cor.mean(S)

shrinkage.intensity(returns, prior, sample)

shrinkage.p(returns, sample)

shrinkage.r(returns, sample, pi.est)

shrinkage.c(prior, sample)

cov.shrink(...)

cov.sample(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{returns}{ A zoo object of returns. This is TxM }
  \item{sample}{ The sample covariance matrix (synonomous to S) }
  \item{prior}{ The shrinkage target covariance matrix (synonomous to F) }
  \item{S}{ The sample covariance matrix }
  \item{pi.est}{ The estimate returned from shrinkage.p }
  \item{\dots}{ Additional parameters to pass to prior.fun }
}
\details{
  cov_shrink(h, prior.fun = cov.prior.cc, ...)

  cov_shrink(h, T, constant.fun, prior.fun = cov.prior.cc, ...)

  cov_shrink(h, ...)


  cov_sample{...}

  cov_sample{returns}

  cov_sample{returns}

  T - Length of returns series used in scaling of shrinkage coefficient

  h - A generic tawny object representing either a returns, covariance, or correlation matrix

  constant.fun - Use this function to calculate the shrinkage constant

  prior.fun - Generates the prior/model covariance matrix 

  Most of the code related to the shrinkage estimator is tied to calculating a
  value for the shrinkage coefficient. The remainder of the code shrinks the
  sample covariance matrix towards the target. In addition, there is a function
  generator used in conjunction with the optimizePortfolio process to produce
  a correlation matrix based on the shrinkage.
}
\value{
  Scalars are produced by all of the shrinkage.* functions, resulting in the 
  final shrinkage coefficient, calculated by shrinkage.intensity.

  The cov.sample function calculates the sample covariance matrix and is MxM.

  The cov.shrink function produces the shrunk version of the covariance matrix
  and has the same dimensions as the sample covariance matrix.

  The cor.mean function calculates the constant correlation used in estimating
  the global mean (aka the shrinkage target) produced by cov.prior.cc.
}
\author{ Brian Lee Yung Rowe }
\seealso{ \code{\link{tawny}}, \code{\link{optimizePortfolio}} }
\examples{
  # Estimate the covariance matrix based on the given asset returns
  data(sp500.subset)
  ys <- TawnyPortfolio(sp500.subset, 150)
  S.hat <- cov_shrink(ys)

  # Optimize the portfolio weights using the shrinkage estimator
  ws <- optimizePortfolio(ys, ShrinkageDenoiser())
  #plotPerformance(ys,ws, bg='white', name='Shrinkage')

  # Calculate the sample covariance matrix
  #S <- cov.sample(ys)

  # Calculate the shrinkage coefficient
  #F <- cov.prior.cc(S)
  #k <- shrinkage.intensity(ys, F, S)
}
\keyword{ ts }
