\name{filter.RMT}
\alias{filter.RMT}
\alias{getRandomMatrix}
\alias{mp.density.hist}
\alias{mp.density.kernel}
\alias{mp.fit.hist}
\alias{mp.fit.kernel}
\alias{mp.eigen.max}
\alias{mp.eigen.min}
\alias{mp.lambdas}
\alias{mp.rho}
\alias{mp.theory}
\alias{cor.empirical}
\alias{r.normalize}
\alias{clean.bouchaud}
\alias{denoise}
\title{ Filter noise from a correlation matrix using RMT to identify the noise }
\description{
  Used to filter out all eigenvalues below k*. At a later date this
  will become pluggable so other people can use their own functions and/or
  provide their own parameters to this function.
}
\usage{
filter.RMT(h, hint, ..., type = 'kernel')

getRandomMatrix(m, t)

mp.theory(Q, sigma, e.values = NULL, steps = 200)

mp.density.hist(h, breaks = NULL, cutoff = 0.01)
mp.density.kernel(h, adjust = 0.2, kernel = 'e', ...)

# Fit the appropriate MP curve to the data. This will estimate Q and sigma.
mp.fit.hist(hist)
mp.fit.kernel(hist)

# Marcenko-Pastur theoretical minimum and maximum eigenvalues
mp.eigen.max(Q, sigma)
mp.eigen.min(Q, sigma)

# Generate eigenvalues for theoretical MP distribution
mp.lambdas(Q, sigma, steps)

# This provides the density of the eigenvalues
mp.rho(Q, sigma, e.values)

r.normalize(h)

cor.empirical(h)

denoise(hist, lambda.plus = 1.6, h = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{ A TxM returns matrix }
  \item{hint}{ A hint to the optimizer for fitting a theoretical curve to the empirical distribution }
  \item{type}{ The type of density calculation to use, either kernel or hist }
  \item{m}{ Scalar representing the number of assets in a portfolio }
  \item{t}{ Scalar representing the number of observations }
  \item{breaks}{ The breaks specification for building a histogram }
  \item{cutoff}{ A threshold to eliminate eigenvalues below this value. Used to improve the performance when correlations are too high and negative eigenvalues exist. }
  \item{adjust}{ Option for density to scale bandwidth }
  \item{kernel}{ Option for density to define what kernel estimator to use }
  \item{hist}{ Histogram-like object that contains eigenvalues, eigenvectors, and density information }
  \item{Q}{ Ratio of data points per time series }
  \item{sigma}{ Standard deviation of entries in h }
  \item{e.values}{ Eigenvalues }
  \item{steps}{ Number of steps to use in plotting theoretical MP distribution }
  \item{lambda.plus}{ The upper eigenvalue bound used to cutoff noisy eigenvalues. The default of 1.6 is based on empirical data and is not optimal. }
  \item{\dots}{ Further arguments to pass to mp.density.* function}
}
\details{
The functions described here are the individual components for filtering a
correlation matrix using Random Matrix Theory. This is only of interest to 
those looking to delve deeper into the mechanics of RMT and/or creating custom
behavior from these building blocks.

In the future it is possible that the portfolio will be kept in a zoo
representation for these internal functions such that all dimensions will be
transposed to T x M.
}
\value{
A cleaned version of the eigenvalues is returned by clean.bouchaud. The 
functions mp.eigen.min/max return scalars while mp.lambdas returns a vector and
mp.rho returns either a scalar or a vector.

}
\author{ Brian Lee Yung Rowe }
\seealso{ \code{\link{tawny}}, \code{\link{optimizePortfolio}} }
\examples{
# Fit the appropriate MP curve to the data. This will estimate Q and sigma.
# hist is a histogram object
h <- getRandomMatrix(1000, 6000)
hist <- mp.density.hist(h)
o <- optim(c(6,1), mp.fit.hist(hist))

# Calculate and plot the theoretical density distribution
rho <- mp.theory(6,1)

clean <- filter.RMT(h, c(6,1))
}
\keyword{ ts }
