\name{prepdata}
\alias{prepdata}
\alias{prepaffy}

\title{Compute useful summary statistics from a data object.}

\usage{
prepdata(expressionset, intgroup, do.logtransform)
prepaffy(expressionset, x)
}

\arguments{
  \item{expressionset}{An object of class
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} for one colour
  non Affymetrix data,
  \code{\link[affy:AffyBatch-class]{AffyBatch}} for Affymetrix data,
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}} for two colour arrays, or
  \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}} for Illumina
  bead arrays.}
  \item{intgroup, do.logtransform}{as in
    \code{\link{arrayQualityMetrics}}.}
  \item{x}{A list, typically the result from a prior call to \code{prepdata}.}
}

\description{\code{prepdata} computes summary statistics that are useful
  for all platforms; \code{prepaffy} computes Affymetrix-specific ones. These are
  helper functions used by \code{\link{arrayQualityMetrics}}.}

\value{A list with various derived quantities. In the case of
  \code{prepaffy}, the returned list is \code{x} with the additional
  elements appended.
}

\details{
See the vignette \emph{Working with arrayQualityMetrics report sections}.
}

\author{Audrey Kauffmann, Wolfgang Huber}
