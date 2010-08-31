\name{aqmobj.box-class}
\docType{class}
\alias{aqmobj.box-class}
\alias{aqmobj.box}
\alias{class.aqmobj.box}
\alias{aqmobj.heat-class}
\alias{aqmobj.heat}
\alias{class.aqmobj.heat}
\alias{aqmobj.ma-class}
\alias{aqmobj.ma}
\alias{class.aqmobj.ma}
\alias{aqmobj.spatial-class}
\alias{aqmobj.spatial}
\alias{class.aqmobj.spatial}
\alias{aqmobj.nuse-class}
\alias{aqmobj.nuse}
\alias{class.aqmobj.nuse}
\alias{aqmobj.rle-class}
\alias{aqmobj.rle}
\alias{class.aqmobj.rle}

\title{ Class to contain data generated from the aqm.boxplot,
  aqm.heatmap, aqm.maplot, aqm.spatial, aqm.nuse and aqm.rle functions. }

\description{ Class to contain data generated from \code{aqm.boxplot},
  \code{aqm.heatmap}, \code{aqm.maplot}, \code{aqm.spatial},
  \code{aqm.nuse} and \code{aqm.rle}. }

\section{Slots}{
  \describe{
    \item{\code{plot}:}{ In the case of \code{aqmobj.box} or
  \code{aqmobj.spatial}, an object of class
      \code{\link[lattice:trellis.object]{trellis.object}} if one
      channel arrays and a list of trellis.object if several channels
  arrays. In the case of \code{aqmobj.heat}, \code{aqmobj.ma} ,
  \code{aqmobj.nuse} or \code{aqmobj.rle}, an object of class
  \code{\link[lattice:trellis.object]{trellis.object}}. }

    \item{\code{section}:}{ A character string with a name for the section
      the plot belongs to in the report. }
    
    \item{\code{title}:}{ A character string with the title of the plot
      to be written in the report. }
    
    \item{\code{legend}:}{ A character string with the legend of the plot
      to be written in the report. }
    
    \item{\code{scores}:}{ A numeric for each array corresponding to the
      scores assessed from the plot. }
    
    \item{\code{outliers}:}{ List or integer of the arrays that are
      outliers using the function \code{boxplot.stats} on the scores. }
    
    \item{\code{shape}:}{ A character "square" or "rect" depending on
    the aspect ratio desired in the report. }
    }
}

\details{
See \code{\link{aqm.prepdata}}, \code{\link{aqm.prepaffy}} and the \code{aqm} Vignette.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.boxplot}}, \code{\link{aqm.heatmap}},
  \code{\link{aqm.maplot}}, \code{\link{aqm.spatial}},
  \code{\link{aqm.nuse}}, \code{\link{aqm.rle}}, \code{\link{aqm.plot}} }

\keyword{classes}
