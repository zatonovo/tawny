\name{divergence}
\alias{divergence}
\alias{divergence.stability}
\alias{divergence.kl}
\alias{divergence_lim}
\alias{stability_lim}
\alias{plotDivergenceLimit.kl}
\title{ Measure the divergence and stability between two correlation matrices }
\description{
  The Kullback-Leibler distance function can be used to measure the divergence
  between two correlation matrices. Although originally designed for probability
  density functions, the literature shows how this can be extended to 
  correlation matrices. By using this function, one can determine objectively
  the effectiveness of a particular filtering strategy for correlation matrices.
}
\usage{
divergence(...)
divergence.kl(sigma.1, sigma.2)
divergence_lim(...)
stability_lim(...)
divergence.stability(h, count, window, filter)
plotDivergenceLimit.kl(m, t.range, ..., overlay = FALSE)
}
\arguments{
  \item{h}{ A zoo object representing a portfolio with dimensions T x M }
  \item{count}{ The number of bootstrap observations to create }
  \item{window}{ The number of samples to include in each observation. Defaults to the anylength of h. }
  \item{filter}{ The correlation filter to measure }
  \item{sigma.1}{ The sample correlation matrix }
  \item{sigma.2}{ The model correlation matrix (aka the filtered matrix) }
  \item{m}{ The number of assets }
  \item{t}{ The number of samples (dates) in an observation }
  \item{t.range}{ A range of date samples. This can be a simple interval so long as it matches the number of samples per asset in the measured correlation matrix. }
  \item{overlay}{ Overlay the divergence limit plot on an existing plot. Default is FALSE. }
  \item{\dots}{ Additional parameters to pass to plot or lines }
}
\details{
divergence(h, count, window = NULL, filter, measure = 'information')
divergence_lim(m, t = NULL)
stability_lim(m, t = NULL)

measure - The type of divergence to calculate. Possible choices are information (default) or stability.
}
\value{
  A summary of the results of the divergence calculation including the mean
  divergence and an effective limit based on a random matrix.
}
\author{ Brian Lee Yung Rowe}
\examples{
# This is autorun outside of examples
tawny:::.init()

data(sp500.subset)
h <- sp500.subset

plotDivergenceLimit.kl(100, 80:499, col='green', ylim=c(0,55))

divergence(h, 25, filter=function(x) denoise(x, create(RandomMatrixFilter)))
divergence(h, 25, filter=function(x) denoise(x, create(ShrinkageFilter)))

}
\keyword{ ts }
