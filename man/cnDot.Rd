\name{cnDot-method}
\alias{cnDot}
\alias{cnDot,catNetwork-method}
\alias{cnDot,catNetwork,character-method}
\alias{cnDot,catNetwork,character-method,character-method}
\alias{cnDot,matrix-method}
\alias{cnDot,matrix,character-method}
\alias{cnDot,matrix,character-method,character-method}
\alias{cnDot,list-method}
\alias{cnDot,list,character-method}
\alias{cnDot,list,character-method,character-method}

\title{Network Description File}

\description{The function generates a dot-file, the native storage format for \code{Graphviz} 
software package, that describes the graph structure of a \code{catNetwork} object.}

\usage{
 cnDot(object, file=NULL, format="ps", style=NULL)
}

\arguments{
 \item{object}{a \code{catNetwork}, a list of \code{catNetwork}s or a parent matrix}
 \item{file}{a \code{character}, an optional output file name}
 \item{format}{a \code{character}, an optional output file format, "ps" or "pdf"}
 \item{style}{a \code{list} of triplets, nodes' shape, color and edge-color}
}

\details{
The function generates a \code{dot}-text file as supported by \code{Graphviz} library. 
In order to draw a graph the user needs a \code{dot}-file converter and \code{pdf}/postscript viewer. 
The environment variables \code{R_DOTVIEWER} and \code{R_PDFVIEWER} specify the corresponding executable routines. 

If \code{Graphviz} is installed and the variable \code{R_DOTVIEWER} is set with the full 
path to the \code{dot} executable file (the routine that converts a \code{dot}-text file to a postscript or \code{pdf}), 
a \code{pdf} or postscript file is created depending on the value of the \code{format} parameter.

If the \code{file} variable is not specified, then the function just prints out the resulting string  
which otherwise would be written into a \code{dot} file. 
Next, if a \code{pdf}-viewer is available, the created postscript or \code{pdf} file is shown.

}

\value{A \code{character} or a \code{dot}-file}

\author{N. Balov}

\examples{
  #cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  #cnDot(object=cnet, file="cnet")
}

\seealso{\code{\link{catnet-package}}, \code{\link{cnPlot}} }

\keyword{methods}
\keyword{graphs}
\keyword{aplot}
