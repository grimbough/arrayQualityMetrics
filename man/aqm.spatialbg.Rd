\name{aqm.spatialbg}
\alias{aqm.spatialbg}

\title{ Performs spatial distribution representation of background
  intensities of the arrays from aqmobj.prepdata objects. }

\description{
  \code{aqm.spatialbg} performs representation of the spatial distribution
  of the background intensities on the arrays, outlier detection and formats the
  output for \code{aqm.plot} usage.
}

\usage{ aqm.spatialbg(expressionset, dataprep, scale) }

\arguments{
  \item{expressionset}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{dataprep}{ An object of class \code{\link{aqmobj.prepdata}} }
  \item{scale}{ The spatial distribution can be represented on the rank
  of the intensities or on the logarithm scale. Possible options are
  thus 'Rank' and 'Log'. }
 }

\value{ An object of class \code{\link{aqmobj.spatialbg}}. }

\details{
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.spatialbg}} }

