\name{aqm.heatmap}
\alias{aqm.heatmap}

\title{ Performs dendrogram on aqmobj.prepdata objects. }

\description{
  \code{aqm.heatmap} performs a dendrogram of the distances between
  arrays, outlier detection from it and formats the output for
  \code{aqm.plot} usage.
}

\usage{ aqm.heatmap(expressionset, dataprep, intgroup = "Covariate", ...) }

\arguments{
  \item{expressionset}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{dataprep}{ An object of class \code{\link{aqmobj.prepdata}} }
  \item{intgroup}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{...}{ Any arguments to \code{\link[lattice:levelplot]{levelplot}} }
}

\value{ An object of class \code{\link{aqmobj.heat}}. }

\details{
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.heat}} }

