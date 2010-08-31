\name{aqmobj.rnadeg-class}
\docType{class}
\alias{aqmobj.rnadeg-class}
\alias{aqmobj.rnadeg}
\alias{class.aqmobj.rnadeg}

\title{ Class to contain data generated from aqm.rnadegplot. }

\description{ Class to contain data generated from \code{aqm.rnadegplot}. }

\section{Slots}{
  \describe{
    \item{\code{plot}:}{ A list to be represented calling the
      \code{\link[affy:AffyRNAdeg]{plotAffyRNAdeg}} function. }
        
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
See \code{\link{aqm.prepdata}} and the \code{aqm} Vignette.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.rnadeg}}, \code{\link{aqm.plot}} }

\keyword{classes}
