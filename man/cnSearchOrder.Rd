\name{cnSearchOrder}
\alias{cnSearchOrder}

\title{Network Search for Given Node Order}

\description{The function implements a MLE based algorithm to search for optimal networks complying 
with a given node order. It returns a list of  networks, with complexities up to some maximal value, that best fit the data. }

\usage{
cnSearchOrder(data, perturbations=NULL, 
	maxParentSet=0, parentSizes=NULL, maxComplexity=0, 
	nodeOrder=NULL, 
	nodeCats=NULL, parentsPool=NULL, fixedParents=NULL, edgeProb=NULL,
	echo=FALSE)
}

\arguments{
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary matrix with the dimensions of \code{data}. A value 1 marks that the node in the corresponding sample as perturbed}
\item{maxParentSet}{an \code{integer}, maximal number of parents for all nodes}
\item{parentSizes}{an \code{integer} vector, maximal number of parents per node}
\item{maxComplexity}{an \code{integer}, the maximal network complexity for the search}
\item{nodeOrder}{a \code{vector} specifying a node order; the search is among the networks consistent with this topological order}
\item{nodeCats}{a \code{list} of node categories}
\item{parentsPool}{a \code{list} of parent sets to choose from}
\item{fixedParents}{a \code{list} of parent sets to choose from}
\item{edgeProb}{a square \code{matrix} of length the number of nodes specifying prior edge probabilities}
\item{echo}{a \code{logical} that sets on/off some functional progress and debug information}
}

\details{
The \code{data} can be a matrix of \code{character} categories with rows specifying the node-variables and columns assumed to be independent samples from an unknown network, or 
a \code{data.frame} with columns specifying the nodes and rows being the samples.

The number of node categories are obtained from the sample. If given, the \code{nodeCats} is used as a list of categories. In that case, \code{nodeCats} should include the node categories presented in the data. 

The function returns a list of networks, one for each admissible complexity within the 
specified range. The networks in the list are the Maximum Likelihood estimates in the class of networks having the given topological order of the nodes and complexity.
When \code{maxComplexity} is not given, thus zero, its value is reset to the maximum possible complexity for the given parent set size. When \code{nodeOrder} is not given or \code{NULL}, the order of the nodes in the data is taken, \code{1,2,...}. 

The parameters \code{parentsPool} and \code{fixedParents} allow the user to put some exclusion/inclusion constrains on the possible parenthood of the nodes. They should be given as lists of index vectors, one for each node. 

The rows in \code{edgeProb} correspond to the nodes in the sample. The [i,j]-th element in \code{edgeProb} specifies 
a prior probability for the j-th node to be a parent of the i-th one. 
In calculating the prior probability of a network all edges are assumed independent Bernoulli random variables. 
The elements of \code{edgeProb} are cropped in the range [0,1], such that the zero probabilities effectively exclude the corresponding edges, while the ones force them. 
}

\value{A \code{catNetworkEvaluate} object}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=12, maxParents=3, numCategories=2)
  psamples <- cnSamples(object=cnet, numsamples=100)
  nodeOrder <- sample(1:12)
  nets <- cnSearchOrder(data=psamples, perturbations=NULL, 
		maxParentSet=2, maxComplexity=36, nodeOrder)
  ## next we find the network with complexity of the original one and plot it
  cc <- cnComplexity(object=cnet)
  cnFind(object=nets, complexity=cc)
}

\seealso{\code{\link{cnSearchSA}} }

\keyword{methods}
\keyword{graphs}
