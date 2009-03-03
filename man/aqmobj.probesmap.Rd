\name{aqmobj.probesmap-class}
\docType{class}
\alias{aqmobj.probesmap-class}
\alias{aqmobj.probesmap}
\alias{class.aqmobj.probesmap}

\title{ Class to contain data generated from aqm.probesmap. }

\description{ Class to contain data generated from aqm.probesmap. }

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
    
    \item{\code{shape}:}{ A character "square" or "rect" depending on
    the aspect ratio desired in the report. }
  }
}

\details{
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{aqm.probesmap}, \code{aqm.plot} }

\keyword{classes}
