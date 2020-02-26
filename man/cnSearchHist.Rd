\name{cnSearchHist}
\alias{cnSearchHist}

\title{Parent Histogram Matrix}

\description{Estimation of the parent matrix of nodes from data. 
The frequency of node edges is obtained by fitting networks consistent to randomly generated node orders.}

\usage{
cnSearchHist(data, perturbations=NULL, 
	maxParentSet=1, parentSizes=NULL, maxComplexity=0, 
	nodeCats=NULL, parentsPool=NULL, fixedParents=NULL, 
	score = "BIC", weight="likelihood", 
	maxIter=32, numThreads=2, echo=FALSE)
}

\arguments{
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary \code{matrix} with the dimensions of \code{data}. A value \code{1} designates the corresponding node in the sample as perturbed}
\item{maxParentSet}{an \code{integer}, the maximal number of parents per node}
\item{parentSizes}{an \code{integer} vector, maximal number of parents per node}
\item{maxComplexity}{an \code{integer}, the maximal network complexity for the search}
\item{nodeCats}{a \code{list} of node categories}
\item{parentsPool}{a list of parent sets to choose from}
\item{fixedParents}{a list of parent sets to choose from}
\item{score}{a \code{character}, network selection score such as "AIC" and "BIC"}
\item{weight}{a \code{character}, specifies how the }
\item{maxIter}{an \code{integer}, the number of single order searches to be performed}	
\item{numThreads}{an \code{integer} value, the number of parallel threads}
\item{echo}{a boolean that sets on/off some functional progress and debug information}
}

\details{The function performs \code{niter} calls of \code{cnSearchOrder} for randomly generated node orders (uniformly over the space of all possible node orders), selects networks according to \code{score} and sum their parent matrices weighted by \code{weight}. Three scoring criteria are currently supported: "BIC", "AIC" and maximum complexity for any other value of \code{score}. The \code{weight} can be 
1) "likelihhod", then the parent matrices are multiplied by the network likelihood, 
1) "score", then the parent matrices are multiplied by the exponential of the network score, 
3) any other value of \code{weihgt} uses multiplier 1. In this case 
the entries in the output \code{matrix} show how many times the corresponding parent-child pairs were found. 


The function can runs \code{numThreads} number of parallel threads each processing different order. 
\code{cnSearchHist} function can be useful for empirical estimation of the relationships in some multivariate categorical data.}

\value{A \code{matrix}}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=8, maxParents=3, numCategories=2)
  psamples  <- cnSamples(object=cnet, numsamples=100)
  mhisto <- cnSearchHist(data=psamples, perturbations=NULL, 
		maxParentSet=2, maxComplexity=20)
  mhisto
}

\seealso{\code{\link{cnMatParents}}, \code{\link{cnSearchOrder}} }

\keyword{methods}
\keyword{graphs}
