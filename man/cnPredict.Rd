\name{cnPredict-method}
\alias{cnPredict}
\alias{cnPredict,catNetwork-method}

\title{Prediction}

\description{Predicts the 'not-available' elements in an incomplete sample.}

\usage{
cnPredict(object, data)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{data}{a data \code{matrix} or \code{data.frame}}
}

\details{\code{Data} should be a matrix or data frame of categorical values or indices. If it is a matrix then the rows should represent \code{object}'s nodes; otherwise, the columns represent the nodes. \code{Data}'s values represent \code{object}'s categories either as characters or indices. 
Indices should be integers in the range from 1 to the number of categories of the corresponding node. 
Prediction is made for those nodes that are marked as not-available (\code{NA}) in the data and is based on maximum probability criterion. For each data instance, the nodes are traversed in their topological order in \code{object} and the categorical values with the maximum probability are assigned.
}

\value{An updated sample matrix}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=3)
  ## generate a sample of size 2 and set nodes 8, 9 and 10 as not-available
  psamples <- matrix(as.integer(1+rbinom(10*2, 2, 0.4)), nrow=10)
  psamples[8, ] <- rep(NA, 2)
  psamples[9, ] <- rep(NA, 2)
  psamples[10, ] <- rep(NA, 2)
  ## make show sample rows are named after the network's nodes
  rownames(psamples) <- cnNodes(cnet)
  ## predict the values of nodes 8, 9 and 10
  newsamples <- cnPredict(object=cnet, data=psamples)
}

\keyword{methods}
\keyword{graphs}
