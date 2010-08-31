\name{aqmobj.dens-class}
\docType{class}
\alias{aqmobj.dens-class}
\alias{aqmobj.dens}
\alias{class.aqmobj.dens}
\alias{aqmobj.pca-class}
\alias{aqmobj.pca}
\alias{class.aqmobj.pca}
\alias{aqmobj.spatialbg-class}
\alias{aqmobj.spatialbg}
\alias{class.aqmobj.spatialbg}
\alias{aqmobj.probesmap-class}
\alias{aqmobj.probesmap}
\alias{class.aqmobj.probesmap}
\alias{aqmobj.qcs-class}
\alias{aqmobj.qcs}
\alias{class.aqmobj.qcs}

\title{ Class to contain data generated from aqm.density, aqm.pca,
  aqm.spatialbg, aqm.probesmap, aqm.qcs. }

\description{ Class to contain data generated from \code{aqm.density},
  \code{aqm.pca}, \code{aqm.spatialbg}, \code{aqm.probesmap}, \code{aqm.qcs}. }

\section{Slots}{
  \describe{
    \item{\code{plot}:}{ In the case of \code{aqmobj.dens} or \code{aqmobj.spatialbg}, an object of class
      \code{\link[lattice:trellis.object]{trellis.object}} if one
      channel arrays and a list of \code{\link[lattice:trellis.object]{trellis.object}} if several channels
  arrays. In the case of \code{aqmobj.pca}, \code{aqmobj.probesmap} or \code{aqmobj.qcs}, an object of class
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
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.density}}, \code{\link{aqm.pca}},
  \code{\link{aqm.spatialbg}}, \code{\link{aqm.probesmap}},
  \code{\link{aqm.qcstats}}, \code{\link{aqm.plot}} }

\keyword{classes}
