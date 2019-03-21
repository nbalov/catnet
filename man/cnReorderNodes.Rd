\name{cnReorderNodes-method}
\alias{cnReorderNodes}
\alias{cnReorderNodes,catNetwork,vector-method}

\title{Reorder Network Nodes}

\description{The function rearranges the nodes of a network according to a new order.}

\usage{
cnReorderNodes(object, nodeIndices)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{nodeIndices}{a \code{vector} representing the new node order}
}

\value{A \code{catNetwork} object.}

\details{Node reordering affects the list of node names, parents and probabilities. 
It is a useful operation in cases when comparison of two networks is needed. }

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnMatParents(cnet)
  cnet1 <- cnReorderNodes(object=cnet, nodeIndices=cnOrder(cnet))
  cnNodes(object=cnet1)
  cnMatParents(cnet1)
}

\keyword{methods}
\keyword{graphs}
