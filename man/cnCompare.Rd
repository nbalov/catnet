\name{cnCompare-method}
\alias{cnCompare}
\alias{cnCompare,catNetwork,catNetwork-method}
\alias{cnCompare,catNetwork,matrix-method}
\alias{cnCompare,catNetwork,list-method}
\alias{cnCompare,catNetwork,catNetworkEvaluate-method}

\title{Network Comparison}

\description{Compares two \code{catNetwork} objects by several criteria}

\usage{
cnCompare(object1, object2, extended = TRUE)
}

\arguments{
\item{object1}{a \code{catNetwork} object}
\item{object2}{a \code{catNetwork} object, \code{matrix}, \code{list} of \code{catNetwork}s or \code{catNetworkEvaluate} object}
\item{extended}{a \code{logical} parameter, specifying whether basic but quicker or extended comparison to be performed}
}

\details{Comparison can be performed only between networks with the same sets of nodes. 
The function considers several topology-related comparison metrics. 

First, directed edge comparison is performed and the true positives (TP), the false positive (FP) and 
the false negatives (FN) are reported assuming \code{object1} to be the 'true' network. 

Second, the difference between the binary parent matrices of the two objects is measured as 
the number of positions at which they differ. This is the so called Hamming distance 
and it is coded as \code{hamm}.
Also, when \code{extended} parameter is set to TRUE, 
the difference between the exponents of the parent matrices is calculated, \code{hammexp}.

Third, the node order difference between the two networks is measured as follows. 
Let us call 'order pair' a pair of indices (i,j) 
such that there is a directed path from j-th node to i-th node in the network, 
which sometimes is denoted by j>i. 
The order comparison is done by counting the false positive and false negative order pairs.

The fourth criteria accounts for the so called 'Markov blanket'.
The term 'Markov pair' is used to denote a pair 
of indices which corresponding nodes have a common child. 
In case of extended comparison, the numbers of false positive and false negative Markov pairs are calculated.

The \code{cnCompare} function returns an object with the following slots:  
1) the number of true positive edges \code{TP}; 
2) the number of false positive edges \code{FP}; 
3) the number of false negative edges \code{FN}; 
4) the F-score, which is the harmonic average of the specificity and sensitivity
5) the number of different elements in the corresponding parent matrices \code{hamm}; 
6) the total number of different elements between all powers of the parent matrices \code{hammexp};

Next three numbers identify the difference in the objects' skeletons (undirected graph structure)

7) the number of true positive undirected edges \code{TP}; 
8) the number of false positive undirected edges \code{FP}; 
9) the number of false negative undirected edges \code{FN}; 

10) the number of false positive order pairs \code{order.fp}; 
11) the number of false negative order pairs \code{order.fn}; 
12) the number of false positive Markov pairs \code{markov.fp}; and  
13) the number of false positive Markov pairs \code{markov.fn}.
It is assumed that the first object represents the ground truth with respect to which the comparison is performed. 

If \code{extended} is set off (FALSE) only the edge (TP, FP, FN) and skeleton (TP, FP, FN) numbers 
are reported, otherwise all distance parameters are calculated. 
Turning off the \code{extended} option is recommended for very large networks (e.g. with number of nodes > 500), 
since the calculation of some of the distance metrics involve matrix calculations for which the function is not optimized 
and can be very slow.
}

\value{A \code{catNetworkDistance} if \code{object2} is \code{catNetwork} and \code{catNetworkEvaluate} otherwise.}

\author{N. Balov}

\examples{
  cnet1 <- cnRandomCatnet(numnodes=10, maxParents=2, numCategories=2)
  cnet2 <- cnRandomCatnet(numnodes=10, maxParents=2, numCategories=2)
  dist <- cnCompare(object1=cnet1, object2=cnet2)
  dist
}

\seealso{\code{\link{catNetworkEvaluate-class}}
}

\keyword{methods}
\keyword{graphs}
