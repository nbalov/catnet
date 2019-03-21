\name{cnNumNodes-method}
\alias{cnNumNodes}
\alias{cnNumNodes,catNetwork-method}

\title{Network Size}

\description{Returns the number of nodes of a \code{catNetwork} object.}

\usage{
cnNumNodes(object)
}

\arguments{
\item{object}{a \code{catNetwork}}
}

\value{an \code{integer}}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnNumNodes(object=cnet)
}

\seealso{\code{\link{cnNodes}} }

\keyword{methods}
\keyword{graphs}
