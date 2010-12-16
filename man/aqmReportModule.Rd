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

\title{Compute quality report modules.}

\description{
  These functions produce objects of class \code{aqmReportModule}
  containing the data (plots, text, and
  if applicable, outlier detection) of the various modules of the
  quality report. Given a list of modules, the report is rendered by the
  \code{\link{aqm.writereport}} function.
}

\usage{
aqm.boxplot(x, subsample = 10000)
aqm.density(x)
aqm.heatmap(x)
aqm.pca(x)
aqm.maplot(x)
aqm.meansd(x)
aqm.probesmap(x)

# Affymetrix specific sections
aqm.qcstats(expressionset)
aqm.pmmm(expressionset)
aqm.rnadeg(expressionset, x)
aqm.rle(x, affyproc)
aqm.nuse(x, affyproc)
}

\arguments{
  \item{x}{The object resulting from a call to \code{\link{prepdata}(expressionset)}.}
  \item{expressionset}{An object of class \code{AffyBatch}.}
   \item{affyproc}{The object resulting from a call to
     \code{\link{prepaffy}(expressionset)}.}
  \item{subsample}{For the boxplot computations, if the arrays have more
    features than indicated by this number, then randomly subsample that
    number of them.}
}

\value{
An object of class \code{aqmReportModule}. For a simple example of the
\code{aqm.*} functions, have a look at the source code of the
\code{aqm.pca} function. 
}

\details{
See the vignette \emph{Working with arrayQualityMetrics report sections}.
}

\author{Audrey Kauffmann, Wolfgang Huber}


