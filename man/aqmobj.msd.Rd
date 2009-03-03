\name{aqmobj.msd-class}
\docType{class}
\alias{aqmobj.msd-class}
\alias{aqmobj.msd}
\alias{class.aqmobj.msd}

\title{ Class to contain data generated from aqm.meansd. }

\description{ Class to contain data generated from aqm.meansd. }

\section{Slots}{
  \describe{
    \item{\code{plot}:}{ A matrix to be represented calling the
      \code{\link[vsn:MeanSdPlot]{MeanSdPlot}} function. }
    
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

\seealso{ \code{aqm.meansd}, \code{aqm.plot} }

\keyword{classes}
