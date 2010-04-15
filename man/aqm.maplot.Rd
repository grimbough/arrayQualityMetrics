\name{aqm.maplot}
\alias{aqm.maplot}

\title{ Performs MA-plots on aqmobj.prepdata objects. }

\description{
  \code{aqm.maplot} performs MA-plots, outlier detection from it and
  formats the output for \code{aqm.plot} usage.
}

\usage{ aqm.maplot(dataprep, ...) }

\arguments{
  \item{dataprep}{ An object of class \code{\link{aqmobj.prepdata}}}
  \item{...}{ Any arguments to
  \code{\link[lattice:panel.smoothScatter]{panel.smoothScatter}} }
}

\value{ An object of class \code{\link{aqmobj.ma}}. }

\details{
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.ma}} }
