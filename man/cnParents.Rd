\name{cnParents-method}
\alias{cnParents}
\alias{cnParents,catNetwork,missing-method}
\alias{cnParents,catNetwork,character-method}
\alias{cnParents,catNetwork,vector-method}

\title{Network Parent Structure}

\description{Returns the list of parents of selected nodes of a \code{catNetwork} object. If \code{which} is not specified, the parents of all nodes are listed.}

\usage{
cnParents(object, which)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{which}{a vector of node indices}
}

\value{A list of named nodes.}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnParents(object=cnet)
}

\seealso{\code{\link{cnMatParents}}, \code{\link{cnEdges}} }

\keyword{methods}
\keyword{graphs}
