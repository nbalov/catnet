\name{cnNodeLoglik}
\alias{cnNodeLoglik}
\alias{cnNodeLoglik,catNetwork-method}

\title{Node Log-likelihood}

\description{For a given data sample, the function calculates the log-likelihood of a node 
with respect to a specified parent set.
}

\usage{
cnNodeLoglik(object, node, data, perturbations=NULL)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{node}{an \code{integer} or a \code{list} of \code{integer}s, node indices in the \code{data}}
\item{data}{a \code{matrix} or \code{data.frame} of categories}
\item{perturbations}{an optional perturbation \code{matrix} or \code{data.frame}}
}

\value{a \code{numeric} value}

\author{N. Balov}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  psamples <- cnSamples(object=cnet, numsamples=100)
  cnNodeLoglik(cnet, node=5, data=psamples)
}

\seealso{\code{\link{cnLoglik}} }

\keyword{methods}
\keyword{graphs}
