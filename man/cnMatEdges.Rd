\name{cnMatEdges-method}
\alias{cnMatEdges}
\alias{cnMatEdges,catNetwork-method}

\title{Network Edge Matrix}

\description{Returns a matrix representing the edges of a \code{catNetwork} object.}

\usage{
cnMatEdges(object)
}

\arguments{
\item{object}{a \code{catNetwork} object}
}

\details{The resulting matrix has two columns and the number of edges rows. Edges are given as ordered pairs of the elements of the first and second columns.}

\value{A matrix of characters.}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnMatEdges(object=cnet)
}

\seealso{\code{\link{cnEdges}}, \code{\link{cnMatParents}} }

\keyword{methods}
\keyword{graphs}
