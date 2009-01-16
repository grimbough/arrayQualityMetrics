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
  \item{...}{ Any arguments to \code{\link[lattice:densityplot]{densityplot}} }
 }

\value{ An object of class \code{\link{aqmobj.probesmap}} }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\examples{
  \dontrun{
    library(CCl4)
    data(CCl4)
    dataprep = aqm.prepdata(CCl4, TRUE)
    pm = aqm.probesmap(dataprep, CCl4)
    }
}