\name{cnNodes-method}
\alias{cnNodes}
\alias{cnNodes,catNetwork,missing-method}
\alias{cnNodes,catNetwork,vector-method}

\title{Netwok Nodes}

\description{Returns the list of nodes of a \code{catNetwork} object.}

\usage{
cnNodes(object, which)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{which}{a vector of node indices}
}

\details{Nodes are represented by characters. When a random \code{catNetwork} object is constructed, 
it takes the default node names \code{N#}, 
where # are node indices. The function returns the node names with indices given by parameter \code{which},
 and all node names if \code{which} is not specified. }

\value{a list of \code{character}s, the node names}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnNodes(object=cnet)
}

\seealso{\code{\link{cnNumNodes}} }

\keyword{methods}
\keyword{graphs}
