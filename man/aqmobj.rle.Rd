\name{aqmobj.rle-class}
\docType{class}
\alias{aqmobj.rle-class}
\alias{aqmobj.rle}
\alias{class.aqmobj.rle}

\title{ Class to contain data generated from aqm.rle. }

\description{ Class to contain data generated from aqm.rle. }

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

    \item{\code{scores}:}{ A numeric for each array corresponding to the
      scores assessed from the plot. }
    
    \item{\code{outliers}:}{ List or integer of the arrays that are
      outliers using boxplot.stats on the scores. }
    
    \item{\code{shape}:}{ A character "square" or "rect" depending on
    the aspect ratio desired in the report. }
    }
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{aqm.rle}, \code{aqm.plot} }

\examples{
    data(MLLaffyprep)
    rle = aqm.rle(MLLaffyprep)
    class(rle)
    show(rle)
    aqm.plot(rle)
}
\keyword{classes}
