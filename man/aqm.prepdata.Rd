\name{prepdata}
\alias{prepdata}
\alias{prepaffy}

\title{Compute useful summary statistics from a data object.}

\usage{
prepdata(expressionset, intgroup, do.logtransform)
prepaffy(expressionset)
}

\arguments{
  \item{expressionset}{An object of class
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} for one colour
  non Affymetrix data,
  \code{\link[affy:AffyBatch-class]{AffyBatch}} for Affymetrix data,
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}} for two colour arrays, or
  \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}} for Illumina
  bead arrays.}
  \item{intgroup, do.logtransform}{as in \code{\link{arrayQualityMetrics}}.}
}

\description{\code{prepdata} computes summary statistics that are useful
  for all platforms; \code{prepaffy} Affymetrix-specific ones. This is a
  helper function used by \code{\link{arrayQualityMetrics}}.}

\value{A list with various derived values.}
\author{Audrey Kauffmann.}

\examples{
   ## See vignette
}
