\name{aqm.rnadeg}
\alias{aqm.rnadeg}

\title{ Performs RNA degradation plot on AffyBatch. }

\description{
  \code{aqm.rnadeg} performs RNA degradation on objects of class
  \code{\link[affy:AffyBatch-class]{AffyBatch}} and formats the output for
  \code{aqm.plot} usage.
}

\usage{ aqm.rnadeg(expressionset, ...) }

\arguments{
  \item{expressionset}{ An object of
    class \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{...}{ Any arguments to \code{\link[affy]{AffyRNAdeg}} }
}

\value{ An object of class \code{\link{aqmobj.rnadeg}}. }

\details{
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqmobj.rnadeg}} }
