\name{aqmReportModule}

\alias{aqm.boxplot}
\alias{aqm.density}
\alias{aqm.heatmap}
\alias{aqm.pca}
\alias{aqm.maplot}
\alias{aqm.meansd}
\alias{aqm.probesmap}
\alias{aqm.qcstats}
\alias{aqm.pmmm}
\alias{aqm.rnadeg}
\alias{aqm.rle}
\alias{aqm.nuse}

\alias{class:aqmReportModule}
\alias{aqmReportModule}
\alias{aqmReportModule-class}

\title{Compute quality report sections.}

\description{
  These functions produce objects of class \code{aqmReportModule}
  containing the data (plots, text, and
  if applicable, outlier detection) of the various sections of the
  quality report. The report sections are rendered by the
  \code{\link{aqm.writereport}} function.
}

\usage{
aqm.boxplot(x, subsample = 10000)
aqm.density(x, outliers = c())
aqm.heatmap(x)
aqm.pca(x, outliers)
aqm.maplot(x)
aqm.meansd(x)
aqm.probesmap(x)

# Affymetrix specific sections
aqm.qcstats(expressionset)
aqm.pmmm(expressionset)
aqm.rnadeg(expressionset)
aqm.rle(x, affyproc)
aqm.nuse(x, affyproc)
}

\arguments{
  \item{x}{The object resulting from a call to \code{\link{prepdata}(expressionset)}.}
  \item{outliers}{Integer vector indicating the arrays to be highlighted as 
    outliers in the plots.}
  \item{subsample}{For the boxplot computations, if the arrays have more
    features than indicated by this number, then randomly subsample that
    number of them.}
  \item{expressionset}{An object of class \code{AffyBatch}.}
   \item{affyproc}{The object resulting from a call to
     \code{\link{prepaffy}(expressionset)}.}
}

\value{
An object of class \code{aqmReportModule}. For a simple example of the
\code{aqm.*} functions, have a look at the source code of the
\code{aqm.pca} function. 
}

\details{
See the \code{aqm} Vignette.
}

\examples{
  ## See vignette
}

\author{Audrey Kauffmann, Wolfgang Huber}


