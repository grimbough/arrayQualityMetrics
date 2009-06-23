\name{aqmobj.qcs-class}
\docType{class}
\alias{aqmobj.qcs-class}
\alias{aqmobj.qcs}
\alias{class.aqmobj.qcs}

\title{ Class to contain data generated from aqm.qcs. }

\description{ Class to contain data generated from aqm.qcs. }

\section{Slots}{
  \describe{
    \item{\code{plot}:}{ An object of class
      \code{\link[lattice:trellis.object]{trellis.object}}. }

    \item{\code{section}:}{ A character string with a name for the section
      the plot belongs to in the report. }
    
    \item{\code{title}:}{ A character string with the title of the plot
      to be written in the report. }
    
    \item{\code{legend}:}{ A character string with the legend of the plot
      to be written in the report. }
    
    \item{\code{shape}:}{ A character "square" or "rect" depending on
    the aspect ratio desired in the report. }

    }  
}
\details{
See the aqm.qcstats help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{aqm.qcstats}, \code{aqm.plot} }
\keyword{classes}