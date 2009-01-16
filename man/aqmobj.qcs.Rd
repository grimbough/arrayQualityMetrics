\name{aqmobj.qcs-class}
\docType{class}
\alias{aqmobj.qcs-class}
\alias{aqmobj.qcs}
\alias{class.aqmobj.qcs}

\title{ Class to contain data generated from aqm.qcs. }

\section{Slots}{
  \describe{
    \item{\code{plot}:}{ An object of class
      \code{\link[lattice:trellis.object]{trellis.object}}. }

    \item{\code{type}:}{ A character string with a name for the section
      the plot belongs to in the report. }
    
    \item{\code{title}:}{ A character string with the title of the plot
      to be written in the report. }
    
    \item{\code{legend}:}{ A character string with the legend of the plot
      to be written in the report. }

    }  
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{aqm.qcstats}, \code{aqm.plot} }
\examples{
    library(ALLMLL)
    data(MLL.A)
    qc = aqm.qcstats(MLL.A)
    class(qc)
    show(qc)
    aqm.plot(qc)
}
\keyword{classes}
