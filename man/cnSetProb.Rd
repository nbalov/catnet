\name{cnSetProb-method}
\alias{cnSetProb}
\alias{cnSetProb,catNetwork-method}
\alias{cnSetProb,catSampleNetwork-method}

\title{Set Probability from Data}

\description{The function sets the probability structure of a network from data according to the Maximum Likelihood criterion.}

\usage{
cnSetProb(object, data, perturbations=NULL, nodeCats=NULL)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{data}{a data \code{matrix} or \code{data.frame}}
\item{perturbations}{a binary matrix with the dimensions of \code{data}}
\item{nodeCats}{a \code{list} of node categories}
}

\details{
The function generates a new probability table for \code{object} and returns an updated \code{catNetwork}.
The graph structure of the \code{object} is kept unchanged. 

The \code{data} can be a \code{matrix} in the node-rows format, 
or a \code{data.frame} in the node-column format. 
If given, the \code{nodeCats} is used as a list of categories. In that case, \code{nodeCats} should include the node categories presented in the data. 
}

\value{\code{catNetwork}}

\author{N. Balov}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=3)
  psamples <- matrix(as.integer(1+rbinom(10*100, 2, 0.4)), nrow=10)
  rownames(psamples) <- cnet@nodes
  newcnet <- cnSetProb(object=cnet, data=psamples)
}

\keyword{methods}
\keyword{graphs}
