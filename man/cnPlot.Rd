\name{cnPlot-method}
\alias{cnPlot}
\alias{cnPlot,catNetwork-method}
\alias{cnPlot,catNetworkEvaluate-method}

\title{Plot Network}

\description{Draws the graph structure of \code{catNetwork} object or 
some diagnostic plots associated with a \code{catNetworkEvaluate}}

\usage{
cnPlot(object, file=NULL)
}

\arguments{
\item{object}{\code{catNetwork} or \code{catNetworkEvaluate} object}
\item{file}{a file name}
}

\details{
First we consider the case when \code{object} is a \code{catNetwork}. 
There are two visualization options implemented - one using 'igraph' and the other 'Graphviz'.
The usage of these two alternatives is controlled by two environment variables - 
the logical one \code{R_CATNET_USE_IGRAPH} and the character one \code{R_DOTVIEWER}, correspondingly. 
If \code{igraph} is installed and \code{R_CATNET_USE_IGRAPH} is set to TRUE, 
the function constructs an \code{igraph} compatible object corresponding to the \code{object} and plot it. 

If \code{igraph} is not found, 
the function generates a \code{dot}-file with name \code{file.dot}, if \code{file} is specified, or \code{unknown.dot} otherwise.
Furthermore, provided that \code{Graphviz} library is found and 
\code{R_DOTVIEWER} points to the \code{dot}-file executable, 
the created earlier \code{dot}-file will be compiled to \code{pdf} or postscript, if \code{object} is a list. 
Finally, if the system has \code{pdf} or postscript rendering capabilities and 
\code{R_PDFVIEWER} variable shows the path to the \code{pdf}-rendering application, the resulting 
\code{pdf}-file will be shown.

In case \code{object} is of class \code{catNetworkEvaluate}, then the function draws six relevant plots: 
likelihood vs. complexity, Hamming (\code{hamm}) and exponential Hamming (\code{hammexp}) distances, 
Markov neighbor distance (FP plus FN), 
and the false positive (\code{fp}) and false negative (\code{fn}) edges vs. complexity.
}

\value{A R-plot or \code{dot}-file or \code{pdf}-file.}

\seealso{\code{\link{cnDot}}, \code{\link{catNetworkEvaluate-class}}, \code{\link{cnCompare}}
}

\author{N. Balov}
\examples{
  ## Set R_CATNET_USE_IGRAPH to TRUE if you want to use 'igraph'
  #Sys.setenv(R_CATNET_USE_IGRAPH=FALSE)
  #cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  #cnPlot(object=cnet)
}

\keyword{methods}
\keyword{graphs}
\keyword{aplot}

