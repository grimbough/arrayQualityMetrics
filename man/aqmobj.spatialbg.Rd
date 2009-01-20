\name{aqmobj.spatialbg-class}
\docType{class}
\alias{aqmobj.spatialbg-class}
\alias{aqmobj.spatialbg}
\alias{class.aqmobj.spatialbg}

\title{ Class to contain data generated from aqm.spatialbg. }

\description{ Class to contain data generated from aqm.spatialbg. }

\section{Slots}{
  \describe{
        \item{\code{plot}:}{ An object of class
      \code{\link[lattice:trellis.object]{trellis.object}} if one
      channel arrays and a list of trellis.object if several channels arrays. }

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

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{aqm.spatialbg}, \code{aqm.plot} }

\examples{
\dontrun{
    library(CCl4)
    data(CCl4)
    dataprep = aqm.prepdata(CCl4, TRUE)
    spb = aqm.spatialbg(CCl4, dataprep)
    class(spb)
    show(spb)
    aqm.plot(spb)
}
}
\keyword{classes}
