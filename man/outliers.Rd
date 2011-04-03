\name{outliers}
\alias{outliers}
\alias{boxplotOutliers}

\title{Helper functions for outlier detection and reporting in arrayQualityMetrics}

\description{
  For an overview over outlier detection, please see the
  corresponding section in the vignette \emph{Advanced topics: Customizing arrayQualityMetrics
    reports and programmatic processing of the output}.
  These two functions are helper functions used by the different report
  generating functions, such as \code{\link{aqm.boxplot}}.
 }

\usage{
outliers(exprs, method = c("KS", "sum", "upperquartile"))
boxplotOutliers(x, coef = 1.5) 
}

\arguments{
  \item{exprs}{A matrix whose columns correspond to arrays, rows to the
    array features.}
  \item{method}{A character string specifying the summary statistic to
    be used for each column of \code{exprs}. See Details.}
  \item{x}{A vector of real numbers.}
  \item{coef}{A number is called an outlier if it is larger than the
    upper hinge plus \code{coef} times the interquartile range.
    Upper hinge and interquartile range are computed by
    \code{\link[stats:fivenum]{fivenum}}.}
 }

\value{For \code{outliers}, an object of class \code{\link{outlierDetection}}.
  For \code{boxplotOutliers}, a list with two elements:
  \code{thresh}, the threshold against which \code{x} was compared, and
  \code{outliers}, an integer vector of indices.
}

\details{
  \code{outliers}: with argument \code{method="KS"}, the function first computes for each array
  the value of the \code{\link[stats:ks.test]{ks.test}} test statistic
  between its distribution of intensities and the pooled distribution of
  intensities from all arrays.
  With \code{"sum"} and \code{"upperquartile"}, it computes the sum or
  the 75 percent quantile. Subsequently, it calls \code{boxplotOutliers}
  to identify the outlying arrays.

  \code{boxplotOutliers} uses a criterion similar to that used in
  \code{\link[grDevices:boxplot.stats]{boxplot.stats}} to
  detect outliers in a set of real numbers. The main difference is that
  in  \code{findOutliers}, only the outliers to the right
  (i.e. extraordinarily large values) are detected.
}

\author{Wolfgang Huber}


