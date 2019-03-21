\name{cnLoglik-method}
\alias{cnLoglik}
\alias{cnLoglik,catNetwork-method}

\title{Sample Log-likelihood}

\description{Calculate the log-likelihood of a sample with respect to a given \code{catNetwork} object}

\usage{
cnLoglik(object, data, perturbations=NULL, bysample=FALSE)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{data}{a data \code{matrix} given in the column-sample format, or a \code{data.frame} in the row-sample format}
\item{perturbations}{a binary \code{matrix} with the dimensions of \code{data}. A value \code{1} designates the corresponding node in the sample as perturbed.}
\item{bysample}{a \code{logical}}
}

\details{If \code{bysample} is set to TRUE, the function output is a vector of log-likelihoods of the individual sample records. Otherwise, the total average of the log-likelihood of the sample is reported.
}

\value{A numeric value}

\author{N. Balov}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents = 3, numCategories = 2)
  psamples <- cnSamples(object=cnet, numsamples=100)
  cnLoglik(object=cnet, data=psamples)
}

\seealso{\code{\link{cnNodeLoglik}} }

\keyword{methods}
\keyword{graphs}
\keyword{distribution}
