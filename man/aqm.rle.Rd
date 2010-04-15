\name{aqm.rle}
\alias{aqm.rle}

\title{ Performs RLE plot on aqmobj.prepaffy objects. }

\description{
  \code{aqm.rle} performs RLE boxplots and outlier detection from it and
  formats the output for \code{aqm.plot} usage.  
}

\usage{ aqm.rle(affyproc, ...) }

\arguments{
  \item{affyproc}{ An object of class \code{\link{aqmobj.prepaffy}} }
  \item{...}{ Any arguments to \code{\link[affyPLM:PLMset-class]{Mbox}} }
}

\value{ An object of class  \code{\link{aqmobj.rle}} }

\details{
See the aqm.prepaffy help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.prepaffy}}, \code{\link{aqmobj.prepaffy}},
  \code{\link{aqmobj.rle}} }
