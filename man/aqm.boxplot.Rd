\name{aqm.boxplot}
\alias{aqm.boxplot}
\alias{aqm.density}
\alias{aqm.heatmap}
\alias{aqm.pca}

\title{Compute boxplots, density plots, heatmap, PCA biplot on a aqmobj.prepdata object.}

\description{
  The functions do boxplots, density plots, heatmap, PCA biplot and
  produce output objects ready for rendering by \code{aqm.plot}.
}

\usage{
aqm.boxplot(expressionset, dataprep, intgroup, subsample = 10000, ...)
aqm.density(expressionset, dataprep, intgroup, outliers, ...)
aqm.heatmap(expressionset, dataprep, intgroup, ...)
aqm.pca(expressionset, dataprep, intgroup, outliers, ...)
}
\arguments{
  \item{expressionset, dataprep, intgroup}{Same as for the function
    \code{\link{arrayQualityMetrics}} }
  \item{outliers}{Indicates arrays to be highlighted as potential
    outliers in the plots.}
  \item{subsample}{For the boxplot computations, if the arrays have more
    features than indicates by this number, then randomly subsample that
    number of them. }
  \item{...}{These arguments get passed on to
    \code{\link[lattice:xyplot]{bwplot}} (in the case of
    \code{aqm.boxplot}),
    \code{\link[lattice:xyplot]{xyplot}} (in the case of
    \code{aqm.density}),
    \code{\link[lattice:xyplot]{levelplot}} (in the case of
    \code{aqm.heatmap})
  }
}

\value{
An object of class \code{\link{aqmobj.box}},
  \code{\link{aqmobj.dens}} or \code{\link{aqmobj.heat}}, respectively.
}

\details{
See \code{\link{aqm.prepdata}} and the \code{aqm} Vignette.
}

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.prepdata}}, \code{\link{aqmobj.prepdata}},
  \code{\link{aqmobj.box}} }

