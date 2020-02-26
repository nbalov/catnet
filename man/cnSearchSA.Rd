\name{cnSearchSA}
\alias{cnSearchSA}

\title{Stochastic Network Search}

\description{This function provides a MLE based network search in the space of node orders by Simulated Annealing. For a given sample from an unknown categorical network, it returns a list of \code{catNetwork} objects, with complexity up to some maximal value, that best fit the data.}

\usage{
cnSearchSA(data, perturbations, 
	maxParentSet=0, parentSizes=NULL, 
	maxComplexity=0, nodeCats=NULL, 
	parentsPool=NULL, fixedParents=NULL, 
	edgeProb=NULL, dirProb=NULL, 
	selectMode = "BIC", 
	tempStart=1, tempCoolFact=0.9, tempCheckOrders=10, maxIter=100, 
	orderShuffles=1, stopDiff=0, 
	numThreads=2, priorSearch=NULL, echo=FALSE)
}

\arguments{
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary matrix with the dimensions of \code{data}. A value \code{1} designates the node in the corresponding sample as perturbed}
\item{maxParentSet}{an \code{integer}, maximal number of parents for all nodes}
\item{parentSizes}{an \code{integer} vector, maximal number of parents per node}
\item{maxComplexity}{an \code{integer}, maximal network complexity for the search}
\item{nodeCats}{a \code{list} of node categories}
\item{parentsPool}{a \code{list} of parent sets to choose from}
\item{fixedParents}{a \code{list} of fixed parent sets}
\item{edgeProb}{a square \code{matrix} of length the number of nodes specifying prior edge probabilities}
\item{dirProb}{a square \code{matrix} of length the number of nodes specifying prior directional probabilities}
\item{selectMode}{a \code{character}, optimization network selection criterion such as "AIC" and "BIC"}
\item{tempStart}{a \code{numeric} value, the initial temperature for the annealing}	
\item{tempCoolFact}{a \code{numeric} value, the temperature multiplicative decreasing factor}
\item{tempCheckOrders}{an \code{integer}, the number of iteration, orders to be searched, with constant temperature}
\item{maxIter}{an \code{integer}, the total number of iterations, thus orders, to be processed}
\item{orderShuffles}{a \code{numeric}, the number of shuffles for generating new candidate orders from the last accepted}
\item{stopDiff}{a \code{numeric} value, stopping epsilon criterion}
\item{numThreads}{an \code{integer} value, the number of parallel threads}
\item{priorSearch}{a \code{catNetworkEvaluate} object from a previous search}
\item{echo}{a \code{logical} that sets on/off some functional progress and debug information}
}

\details{
The function implements a Simulated Annealing version of the Metropolis algorithm by constructing a Markov chain in the space of node orders. Given a currently selected order, the algorithm tries to improve its likelihood score by exploring its neighborhood. The order score is defined as the likelihood of the selected according to \code{selectMode} network from the set of estimated networks compatible with that order. 

The \code{data} can be a matrix of \code{character} categories with rows specifying the node-variables and columns assumed to be independent samples from an unknown network, or 
a \code{data.frame} with columns specifying the nodes and rows being the samples.

The number of categories for each node is obtained from the data. It is the user responsibility to make sure the data can be categorized reasonably. If the data is numerical it will be forcibly coerced to integer one, which however may result to NA entries or too many node categories per some nodes, and in either case to the function failure. Use \code{cnDiscretize} to convert numeric data into categorical. 
If given, the \code{nodeCats} is used as a list of categories. In that case, \code{nodeCats} should include the node categories presented in the data. 

The function returns a list of networks, one for any possible complexity within the specified range. Stochastic optimization, based on the criterion of maximizing the likelihood, is carried on the network with complexity closest to, but not above, \code{maxComplexity}. 
If \code{maxComplexity} is not specified, thus the function is called with the default zero value, then \code{maxComplexity} is set to be the complexity of a network with all nodes having the maximum, \code{maxParentSet}, the number of parents.
The \code{selectMode} parameter sets the selection criterion for the network upon which the 
maximum likelihood optimization is carried on. "BIC" is the default choice, while any value 
different from "AIC" and "BIC" results in the maximum complexity criterion to be used, 
the one which selects the network with complexity given by \code{maxComplexity}. 

The parameters \code{tempStart}, \code{tempCoolFact} and \code{tempCheckOrders} control the Simulated Annealing schedule.

\code{tempStart} is the starting temperature of the annealing process.

\code{tempCoolFact} is the cooling factor from one temperature step to another. 
It is a number between 0 and 1, inclusively; 
For example, if \code{tempStart} is the temperature in the first step, 
\code{tempStart*tempCoolFact} will be temperature in the second.

\code{tempCheckOrders} is the number of proposals, that is, the candidate orders from the current order neighborhood, to be checked before decreasing the temperature. 
If for example \code{maxIter} is 40 and \code{tempCheckOrders} is 4, 
then 10 temperature decreasing steps will be eventually performed. 

The \code{orderShuffles} parameters controls the extend of the current order neighborhood. A value of zero indicates that random orders should be used as proposals. For positive \code{orderShuffles}'s, a candidate order is obtained from the current one by performing \code{orderShuffles} number of times the following operation: a random position is picked up at random (uniformly) and it is exchanged with the position right up next to it. If \code{orderShuffles} is negative, then the operation is: two positions are picked up at random and their values are exchanged. 

\code{maxIter} is the maximum length of the Markov chain. 

\code{orderShuffles} is a number that controls the extent of the order neighborhoods.
Each new proposed order is obtained from the last accepted one by 
\code{orderShuffles} switches of two node indices. 

\code{stopDiff} is a stopping criterion. If at a current temperature, 
after \code{tempCheckOrders} orders being checked, no likelihood improvement of level at least \code{stopDiff} is found, then the SA stops and the function exists. 
Setting this parameter to zero guarantees exhausting all of the maximum allowed 
\code{maxIter} order searches.

The function speeds up the Markov Chain by implementing a pre-computing buffer. It runs \code{numThreads} number of parallel threads each of which process a proposed order. If we have more than one acceptance in the batch, the first one is taken as next order selection. The performance boost is more apparent when the Markov chain has a low acceptance rate, in which case the chain can run up to \code{numThreads}-times faster. 

\code{priorSearch} is a result from previous search. This parameters allows a new search to be initiated  
from the best order found so far. Thus a chain of searches can be constructed with varying parameters 
providing greater adaptability and user control.

See the vignettes for more details on the algorithm.}

\value{A \code{catNetworkEvaluate} object.}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=6, maxParents=2, numCategories=2)
  psamples  <- cnSamples(object=cnet, numsamples=100)
  nets <- cnSearchSA(data=psamples, perturbations=NULL, 
		maxParentSet=1, maxComplexity=16)
  cc <- cnComplexity(object=cnet)
  cnFind(object=nets, complexity=cc)
}

\seealso{\code{\link{cnSearchOrder}} }

\keyword{methods}
\keyword{graphs}
