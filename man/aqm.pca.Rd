\name{aqm.pca}
\alias{aqm.pca}

\title{ Performs Principal Component Analysis on aqmobj.prepdata objects. }

\description{
  \code{aqm.pca} performs a PCA of the arrays and formats the output for
  \code{aqm.plot} usage.
}

\usage{ aqm.pca(expressionset, dataprep, intgroup = "Covariate", ...) }

\arguments{
  \item{expressionset}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{dataprep}{ An object of class \code{\link{aqmobj.prepdata}} }
  \item{intgroup}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{...}{ Any arguments to \code{\link[lattice:levelplot]{levelplot}} }
}

\value{ An object of class \code{\link{aqmobj.pca}}. }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.pca}} }

\examples{
   }
