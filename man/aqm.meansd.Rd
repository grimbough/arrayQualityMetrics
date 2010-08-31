\name{aqm.meansd}
\alias{aqm.meansd}

\title{ Performs Mean/SD plot on aqmobj.prepdata objects. }

\description{
  \code{aqm.meansd} performs Mean/SD plot, and
  formats the output for \code{aqm.plot} usage.
}

\usage{ aqm.meansd(dataprep, ...) }

\arguments{
  \item{dataprep}{ An object of class \code{\link{aqmobj.prepdata}} }
  \item{...}{ Any arguments to \code{\link[vsn:meanSdPlot]{meanSdPlot}} }
 }

\value{ An object of class \code{\link{aqmobj.msd}}. }

\details{
See \code{\link{aqm.prepdata}} and the \code{aqm} Vignette.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.msd}} }

