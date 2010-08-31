\name{aqm.rle}
\alias{aqm.rle}
\alias{aqm.nuse}

\title{ Performs RLE and NUSE plots on aqmobj.prepaffy objects. }

\description{
  The functions do RLE boxplots and NUSE boxplots and
  produce output objects ready for rendering by \code{aqm.plot}.
}

\usage{
 aqm.rle(expressionset, dataprep, affyproc, intgroup, subsample = 10000, ...)
 aqm.nuse(expressionset, dataprep, affyproc, intgroup, subsample = 10000, ...) 
 }

\arguments{
  \item{expressionset, dataprep, intgroup}{ Same as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{affyproc}{ An object of class \code{\link{aqmobj.prepaffy}} }
  \item{subsample}{ Same as for the function
  \code{\link{aqm.boxplot}} }
  \item{...}{These arguments get passed on to
    \code{\link[lattice:bwplot]{bwplot}}
  }
}

\value{ An object of class  \code{\link{aqmobj.rle}} or \code{\link{aqmobj.nuse}} }

\details{
See \code{\link{aqm.prepaffy}} and the \code{aqm} Vignette.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.prepaffy}}, \code{\link{aqmobj.prepaffy}},
  \code{\link{aqmobj.rle}}, \code{\link{aqmobj.nuse}} }
