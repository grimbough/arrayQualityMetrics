\name{aqm.rnadeg}
\alias{aqm.rnadeg}

\title{ Computes RNA degradation plot on AffyBatch. }

\description{
  \code{aqm.rnadeg} performs RNA degradation on objects of class
  \code{\link[affy:AffyBatch-class]{AffyBatch}} and formats the output for
  \code{aqm.plot} usage.
}

\usage{ aqm.rnadeg(expressionset) }

\arguments{
  \item{expressionset}{ An object of
    class \code{\link[affy:AffyBatch-class]{AffyBatch}}}
}

\value{ An object of class \code{\link{aqmobj.rnadeg}}. }

\details{
See \code{\link{aqm.prepdata}} and the \code{aqm} Vignette.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqmobj.rnadeg}} }
