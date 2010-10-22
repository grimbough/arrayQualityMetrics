\name{aqm.boxplot}

\alias{aqm.boxplot}
\alias{aqm.density}
\alias{aqm.heatmap}
\alias{aqm.pca}
\alias{aqm.maplot}
\alias{aqm.meansd}
\alias{aqm.probesmap}

\title{Compute plots.}

\description{
  These functions produce objects containing the various plots used in
  the report, which can then be rendered by \code{aqm.plot}.
}

\usage{
aqm.boxplot(x, subsample = 10000)
aqm.density(x, outliers)
aqm.heatmap(x)
aqm.pca(x, outliers)
aqm.maplot(x)
aqm.meansd(x)
aqm.probesmap(x)
}

\arguments{
  \item{x}{An object resulting from a call to \code{\link{prepdata}}.}
  \item{outliers}{Integer vector indicating the arrays to be highlighted as 
    outliers in the plots.}
  \item{subsample}{For the boxplot computations, if the arrays have more
    features than indicated by this number, then randomly subsample that
    number of them.}
}

\value{
An object of class \code{aqmTrellis} or \code{aqmPlotfun}.
}

\details{
See \code{\link{prepdata}} and the \code{aqm} Vignette.
}

\author{Audrey Kauffmann}


