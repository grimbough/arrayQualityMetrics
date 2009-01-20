\name{aqmobj.ma-class}
\docType{class}
\alias{aqmobj.ma-class}
\alias{aqmobj.ma}
\alias{class.aqmobj.ma}

\title{ Class to contain data generated from aqm.maplot. }

\description{ Class to contain data generated from aqm.maplot. }

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

\seealso{ \code{aqm.maplot}, \code{aqm.plot} }

\examples{
    library(CCl4)
    data(CCl4)
    dataprep = aqm.prepdata(CCl4, TRUE)
    ma = aqm.maplot(dataprep)
    class(ma)
    show(ma)
    aqm.plot(ma)
}
\keyword{classes}
