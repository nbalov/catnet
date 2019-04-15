\name{cnComplexity-method}
\alias{cnComplexity}
\alias{cnComplexity,catNetwork-method}
\alias{cnComplexity,catNetwork,integer-method}
\alias{cnKLComplexity}
\alias{cnKLComplexity,catNetwork-method}
\alias{cnKLComplexity,catNetwork-method,integer-method}

\title{Network Complexity}

\description{Returns the complexity of a network}

\usage{
cnComplexity(object, node=NULL, include.unif=TRUE)
cnKLComplexity(object, node=NULL)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{node}{an \code{integer}, node index}
\item{include.unif}{a \code{logical}}
}

\details{Complexity is a network characteristics that depends both on its graphical structure and 
the categorization of its nodes. 

If \code{node} is specified, then the function returns that node complexity, 
otherwise the total complexity of \code{object}, which is the sum of its node complexities, is reported. 
A node complexity is determined by the number of its parents and their categories. 
For example, a node without parents has complexity \code{1}. 
A node with \code{k} parents with respected number of categories \code{c1,c2, ...,ck}, 
has complexity \code{c1*c2*...*ck}. 
Complexity is always a number that is equal or greater than the number of nodes in the network. 
For a network with specified graph structure, its complexity determines the number of parameters needed 
to define its probability distribution and hence the importance of complexity as network characteristic. 

If \code{include.unif} is set to \code{FALSE}
}

\value{An \code{integer}}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnComplexity(object=cnet)
}

\keyword{methods}
\keyword{graphs}
