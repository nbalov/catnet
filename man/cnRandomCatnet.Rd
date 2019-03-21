\name{cnRandomCatnet}
\alias{cnRandomCatnet}

\title{Random Network}

\description{Creates a random \code{catNetwork} with specified number of nodes, 
number of parents and categories per node.}

\usage{
cnRandomCatnet(numnodes, maxParents, numCategories, p.delta1=0.01, p.delta2=0.01)
}

\arguments{
\item{numnodes}{an \code{integer}, the number of nodes}
\item{maxParents}{an \code{integer}, the maximum number of parents per node}
\item{numCategories}{an \code{integer}, the number of categories for each node. It is the function limitation to support only constant number of node categories.}
\item{p.delta1}{a \code{numeric}}
\item{p.delta2}{a \code{numeric}}
}

\details{
A random set of parents, no more than \code{maxParents}, is assigned to each node along with a random conditional probability distribution with values in the union of [p.delta1, 0.5-p.delta2] and [0.5+p.delta2, 1-p.delta1].
Also, each node is assigned a fixed, thus equal, number of categories, \code{numCategories}. 

The function is designed for evaluation and testing purposes only thus lacking much user control over 
the networks it create. Once created with \code{cnRandomCatnet}, a network can be further 
modified manually node by node. However, this requires direct manipulation of the object's slots and may result in 
a wrong network object. It is recommended that after any manual manipulation a call 
is(\code{object}, "catNetwork") is performed to check the object's integrity.
}

\value{A \code{catNetwork} object}

\author{N. Balov}
\examples{
  cnet <- cnRandomCatnet(numnodes=20, maxParents=3, numCategories=2)
}

\seealso{\code{\link{cnNew}} }

\keyword{methods}
\keyword{graphs}
