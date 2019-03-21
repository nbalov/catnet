\name{cnEntropy}
\alias{cnEntropy}
\alias{cnEdgeDistanceKL}
\alias{cnEdgeDistancePearson}
\alias{cnEntropyOrder}

\title{Pairwise Node Entropy}

\description{Calculates the matrix of conditional entropy for each pair of nodes.}

\usage{
cnEntropy(data, perturbations=NULL)
cnEdgeDistanceKL(data, perturbations)
cnEdgeDistancePearson(data, perturbations)
cnEntropyOrder(data, perturbations=NULL)
}

\arguments{
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary \code{matrix} with the dimensions of \code{data}. A value \code{1} designates the corresponding node in the sample as perturbed.}
}

\details{The conditional entropy of node X with respect to Y is defined as -P(X|Y)logP(X|Y), where P(X|Y) is the sample conditional probability, and this is the value at the (X,Y)'th position in the resulting matrix. 
}

\value{A \code{matrix}}

\author{N. Balov}

\seealso{\code{\link{cnParHist}}}

\keyword{methods}
\keyword{graphs}
