\name{catNetworkDistance-class}
\docType{class}
\alias{catNetworkDistance}
\alias{catNetworkDistance-class}

\title{Class "catNetworkDistance"}
\description{This class contains a list of \code{catNetwork}s and it is the output format of \code{cnEvaluate} function}

\section{Slots}{
  \describe{
        \item{\code{hamm}:}{an \code{integer}, the hamming distance between the parent matrices of the found networks and the original network.}
	\item{\code{hammexp}:}{an \code{integer}, the hamming distance between the exponents of the parent matrices.}
	\item{\code{tp}:}{an \code{integer}, the number of true positives directed edges.}
        \item{\code{fp}:}{an \code{integer}, the number of false positives directed edges.}
	\item{\code{fn}:}{an \code{integer}, the number of false negatives directed edges.}
	\item{\code{sp}:}{a \code{numeric}, the specificity.}
	\item{\code{sn}:}{a \code{numeric}, the sensitivity.}
	\item{\code{fscore}:}{a \code{numeric}, the F-score.}
	\item{\code{skel.tp}:}{an \code{integer}, the number of true positives undirected edges.}
        \item{\code{skel.fp}:}{an \code{integer}, the number of false positives undirected edges.}
	\item{\code{skel.fn}:}{an \code{integer}, the number of false negatives undirected edges.}
	\item{\code{order.fp}:}{an \code{integer}, the number of false positive order relations.}
	\item{\code{order.fn}:}{an \code{integer}, the number of false negative order relations.}
	\item{\code{markov.fp}:}{an \code{integer}, the number of false positive Markov pairs.}
	\item{\code{markov.fn}:}{an \code{integer}, the number of false negative Markov pairs.}
	\item{\code{KLdist}:}{a \code{numerical}, the KL distance, currently inactive.}
  }
}

\details{See in the manual of \code{cnCompare} function for description of different distance criteria.
}

\section{Methods}{
  \describe{
	\item{cnPlot}{\code{signature(object="catNetworkDistance")}:
		Draw some distance plots.}
  }
}

\author{N. Balov}

\seealso{\code{\link{catNetwork-class}}, \code{\link{catNetworkEvaluate-class}}, \code{\link{cnCompare}}, \code{\link{cnPlot}}
}

\keyword{classes}
\keyword{graphs}

