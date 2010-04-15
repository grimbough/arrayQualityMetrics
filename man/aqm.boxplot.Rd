\name{aqm.boxplot}
\alias{aqm.boxplot}

\title{ Performs boxplots on aqmobj.prepdata objects. }

\description{
  \code{aqm.boxplot} performs boxplots, outlier detection from it and
  formats the output for \code{aqm.plot} usage.
}

\usage{
aqm.boxplot(expressionset, dataprep, intgroup = "Covariate",
            grouprep = FALSE, ...) } 

\arguments{
  \item{expressionset}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{dataprep}{ An object of class
    \code{\link[arrayQualityMetrics]{aqmobj.prepdata}} }
  \item{intgroup}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{grouprep}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{...}{ Any arguments to \code{\link[lattice:xyplot]{bwplot}} }
}

\value{ An object of class \code{\link{aqmobj.box}}. }

\details{
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.box}} }

