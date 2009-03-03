\name{aqmobj.nuse-class}
\docType{class}
\alias{aqmobj.nuse-class}
\alias{aqmobj.nuse}
\alias{class.aqmobj.nuse}

\title{ Class to contain data generated from aqm.nuse. }

\description{ Class to contain data generated from aqm.nuse. }

\section{Slots}{
  \describe{
    \item{\code{plot}:}{ A matrix to be represented calling the
      \code{\link{aqm.plot}} function. }

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

\details{
See the aqm.prepaffy help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{aqm.nuse}, \code{aqm.plot} }

\keyword{classes}
