\name{aqm.spatial}
\alias{aqm.spatial}

\title{ Performs spatial distribution representation of the arrays from
  aqmobj.prepdata objects. }

\description{
  \code{aqm.spatial} performs representation of the spatial distribution
  of the intensities on the arrays, outlier detection and formats the
  output for \code{aqm.plot} usage.
}

\usage{ aqm.spatial(expressionset, dataprep) }

\arguments{
  \item{expressionset}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{dataprep}{ An object of class \code{\link{aqmobj.prepdata}} }
 }

\value{ An object of class \code{\link{aqmobj.spatial}}. }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.spatial}} }

\examples{
  \dontrun{
    library(CCl4)
    data(CCl4)
    dataprep = aqm.prepdata(CCl4, TRUE)
    sd = aqm.spatial(CCl4, dataprep, TRUE)
    }
}