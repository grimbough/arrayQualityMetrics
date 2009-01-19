\name{aqmobj.box-class}
\docType{class}
\alias{aqmobj.box-class}
\alias{aqmobj.box}
\alias{class.aqmobj.box}

\title{ Class to contain data generated from aqm.boxplot. }

\description{ Class to contain data generated from aqm.boxplot. }

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
    
    \item{\code{scores}:}{ A numeric for each array corresponding to the
      scores assessed from the plot. }
    
    \item{\code{outliers}:}{ List or integer of the arrays that are
      outliers using boxplot.stats on the scores. }
    }
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{aqm.boxplot}, \code{aqm.plot} }

\examples{
    library(CCl4)
    data(CCl4)
    dataprep = aqm.prepdata(CCl4, TRUE)
    bo = aqm.boxplot(dataprep)
    class(bo)
    show(bo)
    aqm.plot(bo)
}

\keyword{classes}
