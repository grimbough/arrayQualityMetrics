\name{outliers}
\alias{outliers}
\alias{ksOutliers}

\title{Functions for determining outliers}

\description{
  \code{outliers} uses a criterion similar to that used in
  \code{\link[grDevices:boxplot.stats]{boxplot.stats}} to
  detect outliers in a set of real numbers. The main difference is that
  in  \code{outliers}, only the outliers to the right
  (i.e. extraordinarily large values) are detected.
  \code{ksOutliers}: see Details.
 }

\usage{
outliers = function (x, coef = 1.5) 
ksOutliers(exprs, subsamp = 300)
}

\arguments{
  \item{x}{A vector of real numbers.}
  \item{coef}{A number is called an outlier if it is larger than the
    upper hinge plus \code{coef} times the interquartile range.
    Upper hinge and interquartile range are computed by
    \code{\link[stats:fivenum]{fivenum}}.
  }
  \item{exprs}{A matrix whose columns correspond to arrays, rows to the array features.}
  \item{subsamp}{If \code{exprs} has more than \code{subsamp} rows, use only
    a random subsample of size \code{subsamp} for calling \code{\link[stats:ks.test]{ks.test}}.}
 }

\value{For \code{outliers}, an integer vector of indices
  (see \code{\link[base:which]{which}}).
For \code{ksOutliers}, a list with two elements: the KS test statistics
\code{s} and the value of \code{outliers(s)}.}

\details{
 \code{ksOutliers} determines outlier arrays based on
  the value of the \code{\link[stats:ks.test]{ks.test}} test statistic between each
  arrays distribution of intensities and the pooled distribution of
 intensities from all arrays.
}

\author{Wolfgang Huber <whuber@embl.de> }


