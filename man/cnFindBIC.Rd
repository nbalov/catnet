\name{cnFindBIC-method}
\alias{cnFindBIC}
\alias{cnFindBIC,list-method}
\alias{cnFindBIC,catNetworkEvaluate-method}

\title{Find Network by BIC}

\description{
	This is a model selection routine that finds a network in a set of networks using the BIC criteria.
}

\usage{
 cnFindBIC(object, numsamples)
}
\arguments{
 \item{object}{A \code{list} of \code{catNetworkNode} objects or \code{catNetworkEvaluate}}
 \item{numsamples}{The number of samples used for estimating \code{object}}
}

\details{The function returns the network with maximal BIC value from a list of networks 
as obtained from one of the search-functions \code{cnSearchOrder}, \code{cnSearchSA} and \code{cnSearchSAcluster}. 
The formula used for the BIC is \code{log(Likelihood) - 0.5*Complexity*log(numNodes)}.}

\value{A \code{catNetwork} object with optimal BIC value.}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=12, maxParents=3, numCategories=2)
  psamples <- cnSamples(object=cnet, numsamples=10)
  nodeOrder <- sample(1:12)
  nets <- cnSearchOrder(data=psamples, perturbations=NULL, 
	maxParentSet=2, maxComplexity=36, nodeOrder)
  bicnet <- cnFindBIC(object=nets, numsamples=dim(psamples)[2])
  bicnet
}

\seealso{\code{\link{cnFindAIC}}, \code{\link{cnFind}} }

\keyword{methods}
\keyword{graphs}
