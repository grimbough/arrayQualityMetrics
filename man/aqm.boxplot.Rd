\name{aqm.boxplot}
\alias{aqm.boxplot}

\title{ Performs boxplots on aqmobj.prepdata objects. }

\description{
  \code{aqm.boxplot} performs boxplots, outlier detection from it and
  formats the output for \code{aqm.plot} usage.
}

\usage{ aqm.boxplot(dataprep, ...) }

\arguments{
  \item{dataprep}{ An object of class
  \code{\link[arrayQualityMetrics:class.aqmobj.prepdata]{aqmobj.prepdata}} }
   \item{...}{ Any arguments to \code{\link[lattice:bwplot]{bwplot}} }
}

\value{ An object of class \code{\link{aqmobj.box}}. }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.box}} }

\examples{
    library(CCl4)
    data(CCl4)
    dataprep = aqm.prepdata(CCl4, TRUE)
    bo = aqm.boxplot(dataprep)
}