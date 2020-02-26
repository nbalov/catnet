\name{cnProb-method}
\alias{cnProb}
\alias{cnPlotProb}
\alias{cnProb,catNetwork-method}
\alias{cnPlotProb,catNetwork-method}

\title{Conditional Probability Structure}

\description{Returns the list of conditional probabilities of nodes specified by \code{which} parameter of a \code{catNetwork} object. Node probabilities are reported in the following format. First, node name and its parents are given, then a list of probability values corresponding to all combination of parent categories (put in brackets) and node categories. 
For example, the conditional probability of a node with two parents, 
such that both the node and its parents have three categories, 
is given by 27 values, one for each of the 3*3*3 combination.
}

\usage{
 cnProb(object, which=NULL)
 cnPlotProb(object, which=NULL)
}

\arguments{
 \item{object}{a \code{catNetwork} object}
 \item{which}{a \code{vector} of indices}
}

\value{A named list of probability tables. 
}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnProb(object=cnet)
  cnPlotProb(object=cnet)
}

\keyword{methods}
\keyword{graphs}
\keyword{distribution}
