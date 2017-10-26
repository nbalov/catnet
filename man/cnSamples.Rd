\name{cnSamples-method}
\alias{cnSamples}
\alias{cnSamples,catNetwork-method}

\title{Samples from Network}

\description{Generates samples from of a \code{catNetwork} object.}

\usage{
cnSamples(object, numsamples = 1, perturbations = NULL, output="frame", 
	as.index=FALSE, naRate=0)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{numsamples}{an \code{integer}, the number of samples to be generated}
\item{perturbations}{a \code{vector}, node perturbations}
\item{output}{a \code{character}, the output format. Can be a \code{data.frame} or \code{matrix}.}
\item{as.index}{a \code{logical}, the output categorical format}
\item{naRate}{a \code{numeric}, the proportion of NAs per sample instance}
}

\details{If the output format is "matrix" then the resulting sample matrix is in row-node format - the rows correspond to the object's nodes while 
the individual samples are represented by columns. 
If the output format is "frame", which is by default, 
the result is a data frame with columns representing the nodes 
and levels the set of categories of the respected nodes.
If \code{as.index} is set to TRUE, the output sample consists of categorical indices, otherwise, and this is by default, of characters specifying the categories. 

A perturbed sample is a sample having nodes with predefined, thus fixed, values. 
Non-perturbed nodes, the nodes which values have to be set, are designated with zeros in the \code{perturbation} vector and their values are generated conditional on the values of their parents.
While the non-zero values in the perturbation vector are carried on unchanged to the output.

If \code{naRate} is positive, then \code{floor(numnodes*naRate)} NA values are randomly placed in each sample instance.
}

\value{A \code{matrix} or \code{data.frame} of node categories as integers or characters}

\author{N. Balov}
\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=3)
  ## generate a sample of size 100 from cnet
  psamples <- cnSamples(object=cnet, numsamples=100, output="frame", as.index=FALSE)
  ## perturbed sample
  nsamples <- 20
  perturbations <- rbinom(10, 2, 0.4)
  ## generate a perturbed sample of size 100 from cnet
  psamples <- cnSamples(object=cnet, numsamples=nsamples, perturbations, as.index=TRUE)
}

\seealso{\code{\link{cnPredict}} }

\keyword{methods}
\keyword{graphs}
