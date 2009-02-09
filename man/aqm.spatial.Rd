\name{aqm.spatial}
\alias{aqm.spatial}

\title{ Performs spatial distribution representation of the arrays from
  aqmobj.prepdata objects. }

\description{
  \code{aqm.spatial} performs representation of the spatial distribution
  of the intensities on the arrays, outlier detection and formats the
  output for \code{aqm.plot} usage.
}

\usage{ aqm.spatial(expressionset, dataprep, scale) }

\arguments{
  \item{expressionset}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{dataprep}{ An object of class \code{\link{aqmobj.prepdata}} }
  \item{scale}{ The spatial distribution can be represented on the rank
  of the intensities or on the logarithm scale. Possible options are
  thus 'Rank' and 'Log'. }
 }

\value{ An object of class \code{\link{aqmobj.spatial}}. }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.spatial}} }

\examples{
    data(CCl4norm)
    data(CCl4prep)
    sd = aqm.spatial(CCl4norm, CCl4prep, scale = "Log")
}