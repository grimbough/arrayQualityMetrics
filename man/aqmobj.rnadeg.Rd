\name{aqmobj.rnadeg-class}
\docType{class}
\alias{aqmobj.rnadeg-class}
\alias{aqmobj.rnadeg}
\alias{class.aqmobj.rnadeg}

\title{ Class to contain data generated from aqm.rnadegplot. }

\section{Slots}{
  \describe{
    \item{\code{plot}:}{ A list to be represented calling the
      \code{\link[affy:plotAffyRNAdeg]{plotAffyRNAdeg}} function. }
        
    \item{\code{type}:}{ A character string with a name for the section
      the plot belongs to in the report. }
    
    \item{\code{title}:}{ A character string with the title of the plot
      to be written in the report. }
    
    \item{\code{legend}:}{ A character string with the legend of the plot
      to be written in the report. }
  }
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{aqm.rnadeg}, \code{aqm.plot} }

\examples{
    library(ALLMLL)
    data(MLL.A)
    rd = aqm.rnadeg(MLL.A)
    class(rd)
    show(rd)
    aqm.plot(rd)
}
\keyword{classes}
