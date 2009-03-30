\name{divergence}
\alias{divergence}
\alias{divergence.kl}
\alias{divergenceLimit.kl}
\alias{plotDivergenceLimit.kl}
\title{ Measure the divergence between two correlation matrices }
\description{
  The Kullback-Leibler distance function can be used to measure the divergence
  between two correlation matrices. Although originally designed for probability
  density functions, the literature shows how this can be extended to 
  correlation matrices. By using this function, one can determine objectively
  the effectiveness of a particular filtering strategy for correlation matrices.
}
\usage{
divergence(h, count, window = NULL, filter = getCorFilter.RMT())
divergence.kl(sigma.1, sigma.2)
divergenceLimit.kl(m, t = NULL)
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
%\details{
%}
\value{
  A summary of the results of the divergence calculation including the mean
  divergence and an effective limit based on a random matrix.
}
\author{ Brian Lee Yung Rowe}
\examples{
data(sp500.subset)

plotDivergenceLimit.kl(100, 80:499, col='green', ylim=c(0,55))

divergence(sp500.subset, 25, filter=getCorFilter.RMT())
divergence(sp500.subset, 25, filter=getCorFilter.Shrinkage())

}
\keyword{ ts }
