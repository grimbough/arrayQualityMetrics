\name{prepdata}
\alias{prepdata}
\alias{prepaffy}

\title{Compute useful summary statistics from a data object.}

\usage{
prepdata(expressionset, intgroup, do.logtransform, usesvg)
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
  \item{intgroup, do.logtransform, usesvg}{as in \code{\link{arrayQualityMetrics}}.}
}

\description{\code{prepdata} computes summary statistics that are useful
  for all platforms; \code{prepaffy} computes Affymetrix-specific ones. These are
  helper functions used by \code{\link{arrayQualityMetrics}}.}

\value{A list with various derived values.}

\details{
See the vignette \emph{Working with arrayQualityMetrics report sections}.
}

\author{Audrey Kauffmann, Wolfgang Huber}
