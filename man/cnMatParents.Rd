\name{cnMatParents-method}
\alias{cnMatParents}
\alias{cnMatParents,catNetwork,missing-method}
\alias{cnMatParents,catNetwork,vector-method}

\title{Network Parent Matrix}

\description{Returns the binary matrix of parent-child relations of a \code{catNetwork} object.}

\usage{
cnMatParents(object, nodeorder)
}

\arguments{
\item{object}{a \code{catNetwork} or \code{catNetworkFit} object}
\item{nodeorder}{an integer \code{vector} specifying the order of the nodes to be taken}
}

\details{The resulting matrix has a value \code{1} at row \code{i} and column \code{j} if 
\code{i}-th node has \code{j}-th node as a parent, and \code{0} otherwise.}
\value{A matrix}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnMatParents(object=cnet)
}

\seealso{\code{\link{cnParents}}, \code{\link{cnMatEdges}} }

\keyword{methods}
