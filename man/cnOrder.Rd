\name{cnOrder-method}
\alias{cnOrder}
\alias{cnOrder,catNetwork-method}
\alias{cnOrder,list-method}

\title{Network Node Order}

\description{The function returns an 
order of the nodes of a network that is compatible with its parent structure.
}

\usage{
cnOrder(object)
}

\arguments{
\item{object}{a \code{catNetwork} or a \code{list} of node parents.}
}

\details{An order is compatible with the parent structure of a network if 
each node has as parents only nodes appearing earlier in that order. 
That such an order exists is guaranteed by the fact that every \code{catNetwork} is a DAG (Directed Acyclic Graph). The result is one order out of, eventually, many possible. }

\value{a list of node indices.}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=20, maxParents=3, numCategories=2)
  cnOrder(object=cnet)
}


\keyword{methods}
\keyword{graphs}
