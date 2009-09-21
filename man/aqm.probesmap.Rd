\name{aqm.probesmap}
\alias{aqm.probesmap}

\title{ Performs probes mapping on aqmobj.prepdata objects. }

\description{
  \code{aqm.probesmap} performs probes mapping, and
  formats the output for \code{aqm.plot} usage.
}

\usage{ aqm.probesmap(expressionset, dataprep, ...) }

\arguments{
  \item{expressionset}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{dataprep}{ An object of class \code{\link{aqmobj.prepdata}} }
  \item{...}{ Any arguments to \code{\link[lattice:histogram]{densityplot}} }
 }

\details{
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\value{ An object of class \code{\link{aqmobj.probesmap}} }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }
