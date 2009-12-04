\name{aqm.density}
\alias{aqm.density}

\title{ Performs density plots on aqmobj.prepdata objects. }

\description{
  \code{aqm.density} performs density curves, outlier detection from it and
  formats the output for \code{aqm.plot} usage. }

\usage{ aqm.density(expressionset, dataprep, intgroup = "Covariate",
grouprep = FALSE, outliers = NULL, ...) }

\arguments{
  \item{expressionset}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{dataprep}{ An object of class
    \code{\link[arrayQualityMetrics]{aqmobj.prepdata}} }
  \item{intgroup}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{grouprep}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{outliers}{ Vector of numbers identifying arrays to be highlighted
 on the plot }
  \item{...}{ Any arguments to \code{\link[lattice:xyplot]{xyplot}} }
}


\value{ An object of class  \code{\link{aqmobj.dens}}. }

\details{
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.dens}} }

