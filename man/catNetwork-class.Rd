\name{catNetwork-class}
\docType{class}
\alias{catNetwork}
\alias{catNetwork-class}

\title{Class "catNetwork"}

\description{This is the base class in the \code{catnet} package for representing Bayesian networks with categorical values.  
It stores both the graph and probability structure of categorical Bayesian networks. 
Technically, \code{catNetwork} is a S4 type of R-class implemented in object-oriented style, with slots representing object components and members for accessing and manipulating class objects. Below we list the slots of \code{catNetwork} and some of its main members along with the functions for creating \code{catNetwork} objects.}

\section{Objects from the Class}{
  Objects can be created by calls of 

  \code{cnNew(nodes, cats, parents, probs)}  

  \code{cnRandomCatnet(numnodes, maxParents, numCategories)}  

  \code{cnCatnetFromEdges(nodes, edges, numCategories)}

  \code{cnCatnetFromSif(file)}

}

\section{Slots}{
  \describe{
	\item{\code{objectName}}{an optional object name of class \code{character}.}
	\item{\code{numnodes}:}{an integer, the number of nodes in the object.}

	\item{\code{nodes}:}{a \code{vector} specifying the node names.}
	\item{\code{parents}:}{a \code{list} specifying the node parents. 
		The list \code{parents} must be the same length as \code{nodes}. 
		Parents are kept as indices in the \code{nodes} vector. 
    	}

	\item{\code{categories}:}{a \code{list} of characters specifying a set of categorical values for each node.
    	}
	
	\item{\code{probabilities}:}{a numerical \code{list} that for each node specifies a discrete probability distribution - 
		the distribution of the node conditional on its parent set. 
		The elements of \code{probabilities} are lists themselves. See \code{cnProb} function for more details.
	    }
	
	\item{\code{maxParents}:}{an integer, the maximum number of node parents.}
	\item{\code{maxCategories}:}{an integer, the maximum number of node categories.}
	\item{\code{meta}:}{an object of class \code{character} storing some meta-data information.}

	\item{\code{nodeComplexity}:}{a numerical \code{vector}, the node complexities.} 
	\item{\code{nodeLikelihood}:}{a numerical \code{vector}, the node likelihoods of the sample being used for estimation.}

	\item{\code{complexity}:}{an \code{integer}, the network complexity} 
	\item{\code{likelihood}:}{a \code{numerical}, the total likelihood of the sample being used for estimation}

	\item{\code{nodeSampleSizes}:}{a numerical \code{vector}, if the object is an estimate, the node sample sizes.} 

  }
}

\section{Methods}{
  \describe{
	\item{cnNew}{\code{signature}(nodes="\code{vector}", 
		cats="\code{list}", parents="\code{list}", probs="\code{list}"):
		Creating a new class object.}

	\item{cnRandomCatnet}{\code{signature}(numnodes="\code{integer}", 
		maxParents="\code{integer}", numCategories="\code{integer}"):
		Creating a random class object.}

	\item{cnCatnetFromEdges}{\code{signature}(nodes="\code{vector}", edges="\code{list}", numCategories="\code{integer}"):
		Deriving a class object from a list of edges.}

	\item{cnCatnetFromSif}{\code{signature}(file="\code{character}"):
		Creating a class object from a file.}

	\item{cnNumNodes}{\code{signature}(object="\code{catNetwork}"):}

	\item{cnNodes}{\code{signature}(object="\code{catNetwork}", which="\code{vector}"):...}

	\item{cnSubNetwork}{\code{signature}(object="\code{catNetwork}", nodeIndices="\code{vector}", 
		indirectEdges="\code{logical}"):...}

	\item{cnReorderNodes}{\code{signature}(object="\code{catNetwork}", nodeIndices="\code{vector}"):...}

	\item{cnParents}{\code{signature}(object="\code{catNetwork}", which="\code{vector}"):...}

	\item{cnMatParents}{\code{signature}(object="\code{catNetwork}", nodeorder="\code{vector}"):...}

	\item{cnEdges}{\code{signature}(object="\code{catNetwork}", which="\code{vector}"):...}

	\item{cnMatEdges}{\code{signature}(object="\code{catNetwork}"):...}

	\item{cnProb}{\code{signature}(object="\code{catNetwork}"):...}

	\item{cnSetProb}{\code{signature}(object="\code{catNetwork}", psamples="\code{matrix}"):...}

	\item{cnPlot}{\code{signature}(object="\code{catNetwork}"):...}

	\item{cnDot}{\code{signature}(object="\code{catNetwork}", file="\code{character}"):...}

	\item{cnSamples}{\code{signature}(object="\code{catNetwork}", nsamples="\code{integer}"):...}

	\item{cnSamplesPert}{\code{signature}(object="\code{catNetwork}", nsamples="\code{integer}", 
		perturbations="\code{matrix}"):...}

	\item{cnOrder}{\code{signature}(object="\code{catNetwork}"):...}

	\item{cnLoglik}{\code{signature}(object="\code{catNetwork}", psamples="\code{matrix}"):...}

	\item{cnComplexity}{\code{signature}(object="\code{catNetwork}"):...}

	\item{cnEvaluate}{\code{signature}(object="\code{catNetwork}", psamples="\code{matrix}", 
		perturbations="\code{matrix}", maxComplexity="\code{integer}"):...}

	\item{cnPredict}{\code{signature}(object="\code{catNetwork}", psamples="\code{matrix}"):...}

	\item{cnCompare}{\code{signature}(object1="\code{catNetwork}", object2="\code{catNetwork}"):...}
  }
}

\details{
   The \code{catNetwork} class provides a comprehensive general structure for representing 
	discrete Bayesian networks by describing both the graph and probability structures. 
Although available for direct access, the class components, its slots, should not be manipulated directly but using the class members instead. A \code{catNetwork} object integrity can always be checked by calling \code{is(object, "catNetwork")}. 
}

\author{N. Balov}

\seealso{\code{\link{cnRandomCatnet}}, \code{\link{cnCatnetFromEdges}}, \code{\link{cnNew}}, 
 	\code{\link{cnNodes}}, \code{\link{cnEdges}}, \code{\link{cnComplexity}}, \code{\link{cnPlot}}
}

\examples{
  set.seed(123)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=2, numCategories=2)
  cnet
}

\keyword{classes}
\keyword{graphs}

