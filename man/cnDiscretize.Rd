\name{cnDiscretize}
\alias{cnDiscretize}

\title{Data Categorization}

\description{Numerical data discretization using empirical quantiles.
}

\usage{
cnDiscretize(data, numCategories, mode="uniform", qlevels=NULL)
}

\arguments{
\item{data}{a numerical \code{matrix} or \code{data.frame}}
\item{numCategories}{an \code{integer}, the number of categories per node}
\item{mode}{a \code{character}, the discretization method to be used, "quantile" or "uniform"}
\item{qlevels}{a list of \code{integer} vectors, the node discretization parameters}
}

\details{
The numerical \code{data} is discretized into given number of categories, \code{numCategories}, using the empirical node quantiles. 
As in all functions of \code{catnet} package that accept data, if the \code{data} parameter is a \code{matrix} then it is organized in the row-node format. If it is a \code{data.frame}, the column-node format is assumed.

The \code{mode} specifies the discretization model.
Currantly, two discretization methods are supported - "quantile" and "uniform", which is the default choice.

The quantile-based discretization method is applied as follows. 
For each node, the sample node distribution is constructed, 
which is then represented by a sum of non-intersecting classes separated by the quantile points of the sample distribution. Each node value is assigned the class index in which it falls into. 

The uniform discretization breaks the range of values of each node into \code{numCategories} equal intervals or of lengths proportional to the corresponding \code{qlevels} values. 

Currently, the function assigns equal number of categories for each node of the data.
}

\value{A \code{matrix} or \code{data.frame} of indices.}

\examples{
  ps <- t(sapply(1:10, function(i) rnorm(20, i, 0.1)))
  dps1 <- cnDiscretize(ps, 3, mode="quantile")
  hist(dps1[1,])
  qlevels <- lapply(1:10, function(i) rep(1, 3))
  qlevels[[1]] <- c(1,2,1)
  dps2 <- cnDiscretize(ps, 3, mode="uniform", qlevels)
  hist(dps2[1,])
}

\seealso{\code{\link{cnSamples}}}

\author{N. Balov, P. Salzman}

