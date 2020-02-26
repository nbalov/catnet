\name{catNetworkEvaluate-class}
\docType{class}
\alias{catNetworkEvaluate}
\alias{catNetworkEvaluate-class}

\title{Class "catNetworkEvaluate"}

\description{This class contains a list of \code{catNetwork}s together with some diagnostic  metrics and information. \code{catNetworkEvaluate} objects are created automatically as result of calling \code{cnEvaluate} or one of the \code{cnSearch} functions. }

\section{Slots}{
  \describe{
	\item{\code{numnodes}:}{an \code{integer}, the number of nodes in the network.}
	\item{\code{numsamples}:}{an \code{integer}, the sample size used for evaluation.}
    	\item{\code{nets}:}{a \code{list} of resultant networks.}
	\item{\code{complexity}}{an \code{integer} vector, the network complexity.} 
	\item{\code{loglik}}{a \code{numerical} vector, the likelihood of the sample being evaluated.}
        \item{\code{hamm}:}{an \code{integer} vector, the hamming distance between the parent matrices of the found networks and the original network.}
	\item{\code{hammexp}:}{an \code{integer} vector, the hamming distance between the exponents of the parent matrices.}
	\item{\code{tp}:}{an \code{integer} vector, the number of true positives directed edges.}
        \item{\code{fp}:}{an \code{integer} vector, the number of false positives directed edges.}
	\item{\code{fn}:}{an \code{integer} vector, the number of false negatives directed edges.}
	\item{\code{sp}:}{a \code{numeric} vector, the specificity.}
	\item{\code{sn}:}{a \code{numeric} vector, the sensitivity.}
	\item{\code{fscore}:}{a \code{numeric} vector, the F-score.}
	\item{\code{skel.tp}:}{an \code{integer} vector, the number of true positives undirected edges.}
        \item{\code{skel.fp}:}{an \code{integer} vector, the number of false positives undirected edges.}
	\item{\code{skel.fn}:}{an \code{integer} vector, the number of false negatives undirected edges.}
	\item{\code{order.fp}:}{an \code{integer} vector, the number of false positive order relations.}
	\item{\code{order.fn}:}{an \code{integer} vector, the number of false negative order relations.}
	\item{\code{markov.fp}:}{an \code{integer} vector, the number of false positive Markov pairs.}
	\item{\code{markov.fn}:}{an \code{integer} vector, the number of false negative Markov pairs.}
	\item{\code{KLdist}:}{a \code{numerical} vector, the KL distance, currently inactive.}
	\item{\code{time}:}{a \code{numerical}, the processing time in seconds.}
  }
}

\details{
The class \code{catNetworkEvaluate} is used to output the result of two functions: \code{cnEvaluate} and \code{cnSearchSA}. 
The usage of it in the first case is explained next.
The complexity and log-likelihood of the networks listed in \code{nets} slots are stored in 
\code{complexity} and \code{loglik} slots. 
Function \code{cnEvaluate} and \code{cnCompare} fills all the slots from \code{hamm} to \code{markov.fn} 
by comparing these networks with a given network.
See in the manual of \code{cnCompare} function for description of different distance criteria.
By calling \code{cnPlot} upon a \code{catNetworkEvaluate} object, some relevant comparison information can be plotted.

When \code{catNetworkEvaluate} is created by calling \code{cnSearchSA} or \code{cnSearchSAcluster} functions,
\code{complexity} and \code{loglik} contains the information not about the networks in the \code{nets} list, 
but about the optimal networks found during the stochastic search process.
Also, the slots from \code{hamm} to \code{markov.fn} are not used.
}

\section{Methods}{
  \describe{
	\item{cnFind}{\code{signature(object="catNetworkEvaluate", complexity="integer")}:
		Finds a network in the list \code{nets} with specific complexity.}
	\item{cnFindAIC}{\code{signature(object="catNetworkEvaluate")}:
		Finds the optimal network according to AIC criterion.}
	\item{cnFindBIC}{\code{signature(object="catNetworkEvaluate")}:
		Finds the optimal network according to BIC criterion.}
	\item{cnPlot}{\code{signature(object="catNetworkEvaluate")}:
		Draw distance plots.}
  }
}

\author{N. Balov}

\seealso{\code{\link{catNetwork-class}}, \code{\link{catNetworkDistance-class}}, \code{\link{cnCompare}}, \code{\link{cnPlot}}
}

\keyword{classes}
\keyword{graphs}

