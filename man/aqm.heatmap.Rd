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

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.heat}} }

\examples{
    data(CCl4norm)
    data(CCl4prep)
    he = aqm.heatmap(CCl4norm, CCl4prep, intgroup = "RNAintegrity")
}