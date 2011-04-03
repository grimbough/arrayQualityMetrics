\name{modulefunctions}

\alias{aqm.boxplot}
\alias{aqm.density}
\alias{aqm.heatmap}
\alias{aqm.pca}
\alias{aqm.maplot}
\alias{aqm.spatial}
\alias{aqm.meansd}
\alias{aqm.probesmap}
% \alias{aqm.qcstats}
\alias{aqm.pmmm}
\alias{aqm.rnadeg}
\alias{aqm.rle}
\alias{aqm.nuse}

\title{Functions for computing quality report modules.}

\description{
  These functions produce objects of class
  \code{\link{aqmReportModule}} representing the various modules of the quality
  report. Given a list of modules, the report
  is then rendered by the \code{\link{aqm.writereport}} function.
}

\usage{
aqm.boxplot(x, subsample=20000, outlierMethod = "KS")
aqm.density(x)
aqm.heatmap(x)
aqm.pca(x)
aqm.maplot(x, subsample=20000, Dthresh=0.15)
aqm.spatial(x, scale="rank", channels = c("M", "R", "G"))
aqm.meansd(x)
aqm.probesmap(x)

# Affymetrix specific sections
% aqm.qcstats(expressionset)
aqm.pmmm(x)
aqm.rnadeg(expressionset, x)
aqm.rle(x, outlierMethod = "KS")
aqm.nuse(x, outlierMethod = "upperquartile")
}

\arguments{
  \item{x}{An object resulting from a call to \code{\link{prepdata}(expressionset)}.}
  \item{expressionset}{An object of class \code{AffyBatch}.}
  \item{subsample}{For efficiency, some computations are performed not
  on the full set of features (which can be hundreds of thousands on
  some arrays), but on a randomly subset whose size is indicated by this
  number.}
  \item{outlierMethod}{As in \code{\link{outliers}}.}
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
An object of class \code{\link{aqmReportModule}}.}

\details{ 
For a simple example of the \code{aqm.*} functions, have a
look at the source code of the \code{aqm.pca} function.  Please see also
the vignette \emph{Advanced topics: Customizing arrayQualityMetrics
reports and programmatic processing of the output}.}

\author{Audrey Kauffmann, Wolfgang Huber}


