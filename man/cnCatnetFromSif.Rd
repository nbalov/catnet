\name{cnCatnetFromSif}

\alias{cnCatnetFromSif}
\alias{cnCatnetFromBif}
\alias{cnCatnetFromSif,character-method}
\alias{cnCatnetFromBif,character-method}

\title{Categorical Network from Simple Interaction File (SIF) and Bayesian Networks Interchange Format (BIF)}

\description{Creates a \code{catNetwork} object from a SIF/BIF file.}

\usage{
cnCatnetFromSif(file, numcats=2)
cnCatnetFromBif(file)
}

\arguments{
\item{file}{a file name}
\item{numcats}{an \code{integer}, the number of node categories}
}

\details{
The function imports a graph structure from a SIF file by assigning equal number \code{numcats} 
of categories for each of its nodes and a random probability model. 
Subsequently, the probability model can be changed by calling \code{cnSetProb} function.
}

\value{A \code{catNetwork} object}

\author{N. Balov}

\seealso{\code{\link{cnNew}}, \code{\link{cnCatnetFromEdges}}, \code{\link{cnSetProb}} }

\keyword{graphs}
