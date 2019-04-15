\name{cnFindAIC-method}
\alias{cnFindAIC}
\alias{cnFindAIC,catNetworkEvaluate-method}
\alias{cnFindAIC,list-method}

\title{Find Network by AIC}

\description{
	This is a model selection routine that finds a network in a set of networks using the AIC criteria.
}

\usage{
 cnFindAIC(object, numsamples)
}
\arguments{
 \item{object}{A \code{list} of \code{catNetwork} objects or \code{catNetworkEvaluate}}
 \item{numsamples}{an \code{integer}}
}

\details{The function returns the network with maximal AIC value from a list of networks 
as obtained from one of the search-functions \code{cnSearchOrder}, \code{cnSearchSA} and \code{cnSearchSAcluster}.
The formula used for the AIC is \code{log(Likelihood) - Complexity}.}

\value{A \code{catNetwork} object with optimal AIC value.}

\author{N. Balov, P. Salzman}

\examples{
library(catnet)
  cnet <- cnRandomCatnet(numnodes=12, maxParents=3, numCategories=2)
  psamples <- cnSamples(object=cnet, numsamples=10)
  nodeOrder <- sample(1:12)
  nets <- cnSearchOrder(data=psamples, perturbations=NULL, 
	maxParentSet=2, maxComplexity=36, nodeOrder)
  aicnet <- cnFindAIC(object=nets)
  aicnet
}

\seealso{\code{\link{cnFind}}, \code{\link{cnFindBIC}} }

\keyword{methods}
\keyword{graphs}
