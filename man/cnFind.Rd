\name{cnFind-method}
\alias{cnFind}
\alias{cnFindKL}
\alias{cnFind,list-method}
\alias{cnFind,catNetworkEvaluate-method}
\alias{cnFindKL,list-method}
\alias{cnFindKL,catNetworkEvaluate-method}

\title{Find Network by Complexity}

\description{
	This is a model selection routine that finds a network in a set of networks for a given complexity.
}

\usage{
 cnFind(object, complexity = 0, alpha=0, factor=1)
 cnFindKL(object, numsamples)
}

\arguments{
 \item{object}{\code{catNetworkEvaluate} or \code{list} of \code{catNetwork}s}
 \item{complexity}{an \code{integer}, target complexity}
 \item{alpha}{a \code{character} or \code{numeric}}
 \item{factor}{a \code{numeric}}
 \item{numsamples}{an \code{integer}}
}

\details{The complexity must be at least the number of nodes of the networks. If no network with the requested complexity exists in the list, then the one with the closest complexity is returned. 
Alternatively, one can apply some standard model selection with \code{alpha}="BIC" and \code{alpha}=AIC. 
}

\value{A \code{catNetwork} object.}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxParents=2, numCategories=2)
  psamples <- cnSamples(object=cnet, numsamples=100)
  netlist <- cnSearchOrder(data=psamples, maxParentSet=2)
  bnet <- cnFind(object=netlist, complexity=cnComplexity(cnet))
  bnet
}

\seealso{\code{\link{cnFindAIC}}, \code{\link{cnFindBIC}} }

\keyword{methods}
\keyword{graphs}
