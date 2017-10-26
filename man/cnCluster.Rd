\name{cnCluster-method}
\alias{cnCluster}
\alias{cnClusterSep}
\alias{cnClusterMI}
\alias{cnCluster,catNetwork-method}
\alias{cnClusterSep,catNetwork-method}

\title{Network Clustering}

\description{
	Retrieving the clusters, the connected sub-networks, of a given network. Estimating the clusters from data.
}

\usage{
 cnCluster(object)
 cnClusterSep(object, data, perturbations=NULL)
 cnClusterMI(data, perturbations=NULL, threshold=0)
}
\arguments{
\item{object}{a \code{catNetwork}}
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary perturbation matrix with the dimensions of \code{data}}
\item{threshold}{a \code{numeric} value}
}

\details{The function \code{cnCluster} constructs a list of subsets of nodes of the \code{object}, each representing a connected sub-network. Isolated nodes, these are nodes not connected to any other, are not reported. Thus, every element of the output list contains at least two nodes. 
The function \code{cnClusterMI} clusters the nodes of the \code{data} using the pairwise mutual information and critical value \code{threshold}.
}

\value{A \code{list} of named nodes.}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=30, maxParents=2, numCategories=2)
  cnCluster(object=cnet)
}

\keyword{methods}
\keyword{graphs}
