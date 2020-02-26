\name{cnEdges-method}
\alias{cnEdges}
\alias{cnEdges,catNetwork,missing-method}
\alias{cnEdges,catNetwork,vector-method}
\alias{cnEdges,catNetwork,character-method}

\title{Network Edges}

\description{Returns the set of directed edges of a \code{catNetwork} object.}

\usage{
cnEdges(object, which)
}
\arguments{
\item{object}{a \code{catNetwork}}
\item{which}{a \code{vector} of node indices or node names}
}

\details{The edges of a \code{catNetwork} are specified as parent-to-child vectors. 
  The function returns a list that for each node with index in the vector \code{which} 
  contains its set of children. 
  If \code{which} is not specified, the children of all nodes are listed.
}

\value{A list of nodes' children.}

\author{N. Balov, P. Salzman}
\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnEdges(object=cnet)
}

\seealso{\code{\link{cnParents}} }

\keyword{methods}
\keyword{graphs}
