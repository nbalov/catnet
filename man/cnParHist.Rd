\name{cnParHist-method}
\alias{cnParHist}
\alias{cnParHist,catNetworkEvaluate-method}
\alias{cnParHist,list-method}

\title{Parenthood Histogram}

\description{Calculates the histogram of parent-child edges for a \code{catNetworkEvaluate} object 
or a list of \code{catNetwork}s}

\usage{
cnParHist(object)
}

\arguments{
\item{object}{\code{catNetworkEvaluate} or list of \code{catNetwork} objects}
}

\value{a numerical \code{matrix}}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=20, maxParents=3, numCategories=2)
  psamples <- cnSamples(cnet, 100)
  nodeOrder <- sample(1:20)
  nets <- cnSearchOrder(psamples, perturbations=NULL, 
	maxParentSet=2, maxComplexity=50, nodeOrder)
  cnParHist(object=nets)
}

\keyword{methods}
\keyword{graphs}
