\name{catnet-package}
\alias{catnet-package}
\alias{catnet}

\docType{package}

\title{catnet: Categorical Bayesian network inference}

\description{
\code{catnet} package provides tools for learning categorical Bayesian networks from data with focus on model selection.
A Bayesian network is defined by a graphical structure in form of directed acyclic graph and a probability model given as a set of conditional distributions, one for each node in the network. 
Considered in the package are only categorical Bayesian networks - networks which nodes represent discrete random variables.
The learning functions implemented in \code{catnet} are based on exhaustive search and output sets of networks with increasing complexity that fit the data according to MLE-based criteria. The final network selection is left to the user. 
These selected networks represent the relations between the node-variables in the data optimally. 

Before starting to use the package, we suggest the user to take a look at some of the main objects used in \code{catnet} such as \code{catNetwork} and \code{catNetworkEvaluate} and then familiarize with the main search functions \code{cnSearchOrder} and \code{cnSearchSA}. 
More details and examples can be found in the manual pages and the vignettes accompanying the package. 

Since \code{catnet} does not have its own plotting abilities, the user needs to setup some external tools in order to visualize networks, or more precisely, \code{catNetwork} objects. 
\code{catnet} provides interface to the \code{Graphviz} library for visualizing graphs. \code{Graphviz} is not a R-package but a platform independent library that the user have to install in advance on its machine in order to use this option.

In order to use \code{Graphviz}, in addition to installing the library, the user has to register a environmental variable with name \code{R_DOTVIEWER} with the path to the \code{Dot} executable file of \code{Graphviz}. The \code{Dot} routine generates a postscript or \code{pdf}-file from a text dot-file. 
Also, the user needs a postscript and \code{pdf}-viewer. The full path to it 
has to be given in another variable with name \code{R_PDFVIEWER}. 
Note that \code{R_PDFVIEWER} variable might be already setup. To check this call \code{Sys.getenv("R_PDFVIEWER")} in R.

The variables \code{R_DOTVIEWER} and eventually \code{R_PDFVIEWER} can be registered in the \code{.First} function residing in the \code{.Rprofile} initializing file.

Below we give two examples. 
On UNIX platform the user may use code like this one

   .First <- function() \{

	......................


  	Sys.setenv(R_DOTVIEWER="/usr/bin/dot")

   \}

On Windows platform the user may have the following two lines in its \code{.First} function 

  .First <- function() \{

	......................

	Sys.setenv(R_PDFVIEWER="\"C:/Program Files (x86)/Adobe/Reader 9.0/Reader/AcroRd32\"")

  	Sys.setenv(R_DOTVIEWER="\"C:/Program Files (x86)/Graphviz 2.26.3/bin/Dot\"")

 \}

  Note that all paths in Windows should be embraced by comment marks, "\"".

}

\author{N. Balov}

