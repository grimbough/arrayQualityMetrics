\name{aqm.density}
\alias{aqm.density}

\title{ Performs density plots on aqmobj.prepdata objects. }

\description{
  \code{aqm.density} performs density curves, outlier detection from it and
  formats the output for \code{aqm.plot} usage.
}

\usage{ aqm.density(dataprep, ...) }

\arguments{
  \item{dataprep}{ An object of class \code{\link{aqmobj.prepdata}} }
  \item{...}{ Any arguments to \code{\link[lattice:xyplot]{xyplot}} }
 }

\value{ An object of class  \code{\link{aqmobj.dens}}. }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.dens}} }

\examples{
    data(CCl4prep)
    de = aqm.density(CCl4prep)
}