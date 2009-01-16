\name{aqmobj.dens-class}
\docType{class}
\alias{aqmobj.dens-class}
\alias{aqmobj.dens}
\alias{class.aqmobj.dens}

\title{ Class to contain data generated from aqm.density. }

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

    }  
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{aqm.density}, \code{aqm.plot} }

\examples{
    library(CCl4)
    data(CCl4)
    dataprep = aqm.prepdata(CCl4, TRUE)
    de = aqm.density(dataprep)
    class(de)
    show(de)
    aqm.plot(de)
}
\keyword{classes}
