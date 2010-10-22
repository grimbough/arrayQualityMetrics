\name{aqm.spatial}
\alias{aqm.spatial}
\alias{aqm.spatialbg}

\title{Computes spatial distribution statistics of the arrays.}

\description{
  \code{aqm.spatial} and \code{aqm.spatialbg} compute spatial distribution statistics of the arrays,
  perform outlier detection and prepare the
  output for \code{aqm.plot} usage.
}

\usage{
  aqm.spatial(expressionset, dataprep, scale)
  aqm.spatialbg(expressionset, dataprep, scale)
}

\arguments{
  \item{expressionset}{ Same input as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{dataprep}{ Result of the call to \code{\link{prepdata}}. }
  \item{scale}{ The spatial distribution can be represented on the rank
  of the intensities or on the logarithm scale. Possible options are
  thus 'Rank' and 'Log'.}
}

\value{
An object of class \code{aqmTrellis}.
}

\details{
See the \code{aqm} vignette.
}

\author{Audrey Kauffmann}

\seealso{ \code{\link{aqm.boxplot}}, \code{\link{prepdata}} }

