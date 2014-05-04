\name{denoise}
\alias{denoise}
\alias{cor.clean}
\alias{cor.empirical}
\alias{normalize}
\alias{deform} % TODO: Move this out
\title{ Remove noise from a correlation matrix using RMT to identify the noise }
\description{
  Used to filter out all eigenvalues below k*. At a later date this
  will become pluggable so other people can use their own functions and/or
  provide their own parameters to this function.
}
\usage{
denoise(...)
cor.empirical(h)
cor.clean(es, lambda.plus=1.6, h=NULL)
normalize(...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{ A returns object }
  \item{es}{ Eigenvalues and vectors }
  \item{lambda.plus}{ Eigenvalue cutoff for cleaning }
  \item{\dots}{ Arguments }
}
\details{
denoise(p, estimator)

normalize(h)

These are different implementations of the denoise function based on the 
estimator provided.
}
\value{
A cleaned correlation matrix.
}
\author{ Brian Lee Yung Rowe }
\seealso{ \code{\link{tawny}}, \code{\link{optimizePortfolio}} }
\examples{
\dontrun{
data(sp500.subset)

h <- sp500.subset
p <- TawnyPortfolio(h, 150)
r1 <- denoise(p, SampleDenoiser())
r2 <- denoise(p, EmpiricalDenoiser())
r3 <- denoise(p, ShrinkageDenoiser())
r4 <- denoise(p, RandomMatrixDenoiser())
}
}
\keyword{ ts }
