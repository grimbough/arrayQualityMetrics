\name{outliers}
\alias{outliers}
\alias{findOutliers}
\alias{class:outlierDetection}

\title{Outlier detection and reporting in arrayQualityMetrics}

\description{
  \code{outliers} uses a criterion similar to that used in
  \code{\link[grDevices:boxplot.stats]{boxplot.stats}} to
  detect outliers in a set of real numbers. The main difference is that
  in  \code{outliers}, only the outliers to the right
  (i.e. extraordinarily large values) are detected.
  \code{ksOutliers}: see Details.
 }

\usage{
findOutliers(x, coef = 1.5) 
outliers(exprs, method = c("KS", "sum", "upperquartile"))
}

\arguments{
  \item{x}{A vector of real numbers.}
  \item{coef}{A number is called an outlier if it is larger than the
    upper hinge plus \code{coef} times the interquartile range.
    Upper hinge and interquartile range are computed by
    \code{\link[stats:fivenum]{fivenum}}.
  }
  \item{exprs}{A matrix whose columns correspond to arrays, rows to the
    array features.}
  \item{method}{A character string specifying the summary statistic to
  be used for each column of \code{exprs}. See details.}
 }

\value{For \code{findOutliers}, a list with two elements:
  \code{thresh}, the threshold against which \code{x} was compared, and
  \code{outliers}, an integer vector of indices.
  For \code{outliers}, an object of class \code{outlierDetection}.}

\details{
 With \code{method="KS"}, the \code{outliers} computes for each array
  the value of the \code{\link[stats:ks.test]{ks.test}} test statistic
  between its distribution of intensities and the pooled distribution of
  intensities from all arrays.
  With \code{"sum"} and \code{"upperquartile"}, it computes the sum or
  the 75 percent quantile.
}

\author{Wolfgang Huber}


