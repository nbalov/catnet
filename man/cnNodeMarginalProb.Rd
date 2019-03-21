\name{cnNodeMarginalProb-method}
\alias{cnNodeMarginalProb}
\alias{cnNodeMarginalProb,catNetwork-method}
\alias{cnJointProb}
\alias{cnJointProb,catNetwork-method}
\alias{cnCondProb}
\alias{cnCondProb,catNetwork-method}

\title{Probability Calculations}

\description{Marginal probability of a node, joint probability of a set of nodes or conditional probability of two sets of nodes.
}

\usage{
cnNodeMarginalProb(object, node)
cnJointProb(object, nodes)
cnCondProb(object, x, y)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{node}{an \code{integer}, a node index in \code{object}}
\item{nodes}{a \code{vector} of node names or indices in \code{object}}
\item{x,y}{\code{vector}s of node categories (either characters or indices) named after nodes of \code{object}}
}

\details{\code{cnJointProb} returns a \code{matrix} with probability values for each combinations of categories arranged in columns. \code{cnCondProb} calculates the value of P(X=x|Y=y).}

\value{a numerical or numerical \code{matrix}}

\author{N. Balov}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnNodeMarginalProb(cnet, node=5)
  cnCondProb(cnet, x=c("N1"=1, "N2"=2), y=c("N3"=1, "N4"=2, "N5"=2))
}

\seealso{\code{\link{cnProb}} }

\keyword{methods}
\keyword{graphs}
