\name{cnCatnetFromEdges}
\alias{cnCatnetFromEdges}
\alias{cnCatnetFromEdges,character-method}

\title{catNetwork from Edges}

\description{Creates a \code{catNetwork} object from list of nodes and edges.}

\usage{
cnCatnetFromEdges(nodes, edges, numCategories=2)
}

\arguments{
\item{nodes}{a \code{vector} of node names}
\item{edges}{a \code{list} of node edges}
\item{numCategories}{an \code{integer}, the number of categories per node}
}

\details{
The function uses a list of nodes and directional edges to create a \code{catNetwork} with specified (fixed) number of node categories.
A random probability model is assigned, which can be changed later by \code{cnSetProb} for example. 
Note that \code{cnSetProb} takes a given data sample and changes both the node categories and their conditional probabilities 
according to it.
}

\value{A \code{catNetwork} object}

\author{N. Balov}

\seealso{\code{\link{cnNew}}, \code{\link{cnCatnetFromSif}}, \code{\link{cnSetProb}}}

\keyword{graphs}
