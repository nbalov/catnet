\name{isDAG}
\alias{isDAG}

\title{Check Direct Acyclic Graph (DAG) Condition}

\description{For a pair of node and parent lists, the function checks whether the DAG condition holds or not.}

\usage{
isDAG(lnodes, lpars)
}

\arguments{
\item{lnodes}{a list of nodes}
\item{lpars}{a list of node parents}
}

\details{The DAG verification algorithm is based on the topological ordering of the graph nodes. If node ordering is not possible, the graph is not a DAG.}

\value{A logical TRUE/FALSE value.}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=20, maxParents=3, numCategories=2)
  isDAG(lnodes=cnet@nodes, lpars=cnet@parents)
}
\keyword{graphs}
