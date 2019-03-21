\name{cnSubNetwork-method}
\alias{cnSubNetwork}
\alias{cnSubNetwork,catNetwork-method}
\alias{cnSubNetwork,catNetwork,vector,logical-method}

\title{Sub-Network}

\description{Returns a sub-network of a given \code{catNetwork} object.}

\usage{
cnSubNetwork(object, nodeIndices, indirectEdges)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{nodeIndices}{a \code{vector}, the subset of nodes to be taken}
\item{indirectEdges}{a \code{logical}, should the indirect connectivity be preserved}
}

\value{A \code{catNetwork} object.}

\details{The function creates a new network from a given one using a subset of its nodes, specified by \code{nodeIndices}. 
If \code{indirectIndices} is set to TRUE, 
then the resulting network contains edges between all nodes 
that are connected by chains of directed edges in the original one. 
The default value of \code{indirectIndices} is FALSE, 
thus the new set of edges is subset of the original one.}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnet1 <- cnSubNetwork(object=cnet, nodeIndices=c(1,2,3,4,5), indirectEdges=TRUE)
  cnNodes(object=cnet)
  cnNodes(object=cnet1)
}

\keyword{methods}
\keyword{graphs}
