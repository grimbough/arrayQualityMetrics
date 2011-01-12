\name{aqmReportModule}

\alias{aqm.boxplot}
\alias{aqm.density}
\alias{aqm.heatmap}
\alias{aqm.pca}
\alias{aqm.maplot}
\alias{aqm.spatial}
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
  These functions produce objects of class
  \code{aqmReportModule} representing the various modules of the quality
  report (see Section \emph{Value}). Given a list of modules, the report
  is rendered by the \code{\link{aqm.writereport}} function.
}

\usage{
aqm.boxplot(x, subsample=10000)
aqm.density(x)
aqm.heatmap(x)
aqm.pca(x)
aqm.maplot(x, subsample=1000, Dthresh=0.15)
aqm.spatial(x, scale="rank", channels = c("M", "R", "G"))
aqm.meansd(x)
aqm.probesmap(x)

# Affymetrix specific sections
aqm.qcstats(expressionset)
aqm.pmmm(x)
aqm.rnadeg(expressionset, x)
aqm.rle(x)
aqm.nuse(x)
}

\arguments{
  \item{x}{An object resulting from a call to \code{\link{prepdata}(expressionset)}.}
  \item{expressionset}{An object of class \code{AffyBatch}.}
  \item{subsample}{For efficiency, some computations are performed not
  on the full set of features (which can be hundreds of thousands on
  some arrays), but on a randomly subset whose size is indicated by this
  number.}
  \item{Dthresh}{In \code{maplot}, the arrays with a Hoeffding D
  statistic larger than this value are called \emph{outliers}. See also
  \code{\link[Hmisc:hoeffd]{hoeffd}}.}
  \item{scale, channels}{In \code{aqm.spatial}, \code{scale} determines the choice of the false colour scale in the
  spatial plots. If the value is \code{"rank"}, then the colour is proportional to the
  ranks of the values; if it is \code{"direct"}, then it is proportional
  to the values themselves. \code{channels} determines for which elements of
  \code{x} spatial plots are made.}
 }


\value{
An object of class \code{aqmReportModule}. The object contains the code
needed to render the plot, a textual legend, a title, the name of
the section that the module should be part of, and, if applicable,
the results of the outlier detection for this modules.
The interpretation of the plot is explained in the legend.}

\details{
For a simple example of the \code{aqm.*} functions, have a look at the
source code of the \code{aqm.pca} function.
Some, but not all of these modules, provide outlier detection. These are:
\code{aqm.heatmap},
\code{aqm.boxplot},
\code{aqm.maplot},
\code{aqm.spatial},
\code{aqm.rle},
\code{aqm.nuse}.
See also the vignette \emph{Customizing arrayQualityMetrics reports}.
}

\author{Audrey Kauffmann, Wolfgang Huber}


