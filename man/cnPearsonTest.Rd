\name{cnPearsonTest-method}
\alias{cnPearsonTest}
\alias{cnPearsonTest,catNetwork-method}

\title{Goodness of Fit Test}

\description{The function calculates the Pearson's chi-square statistics for all nodes of a network.}

\usage{
cnPearsonTest(object, data)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{data}{a data \code{matrix} or \code{data.frame}}
}

\details{For given \code{data} and network \code{object}, the function reports both the chi-square statistics and the degree of freedom for each node in the network for the purpose of performing goodness of fit tests. 
}

\value{A \code{list}}

\author{N. Balov}

\keyword{methods}
\keyword{graphs}
