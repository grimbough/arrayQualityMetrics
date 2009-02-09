\name{aqmobj.prepaffy-class}
\docType{class}
\alias{aqmobj.prepaffy-class}
\alias{aqmobj.prepaffy}
\alias{class.aqmobj.prepaffy}

\title{ Class to contain data generated from aqm.prepaffy. }

\description{ Container for the output of \code{\link{aqm.prepaffy}} and
  for the input of the \code{\link{aqm.rle}} and \code{\link{aqm.nuse}}
  functions. }

\section{Slots}{
  \describe{
    \item{\code{dataPLM}:}{ A \code{\link[affyPLM:PLMset]{PLMset}}. }
    \item{\code{sN}:}{ Integers numbering the arrays to be used to label
      the plots. }
    }
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\examples{
    library(ALLMLL)
    data(MLL.A)
    MLLaffyprep = aqm.prepaffy(MLL.A, sampleNames(MLL.A))
    class(MLLaffyprep)
}
\keyword{classes}
