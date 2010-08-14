\name{ksOutliers}
\alias{ksOutliers}

\title{Function for outlier determing potential outliers}

\description{
  \code{ksOutliers} determines potential outlier arrays based on
  the value of the \code{\link[stats:ks.test]{ks.test}} test statistic between each
  array and the pooled distribution of all arrays.
}

\usage{ksOutliers(x, subsamp = 300, theta = 2)}

\arguments{
  \item{x}{A matrix whose columns correspond to arrays, rows the array features.}
  \item{subsamp}{If \code{x} has more than \code{subsamp} rows, use only
    a random subsample of size \code{subsamp} for calling \code{\link[stats:ks.test]{ks.test}}.}
  \item{theta}{An array is called an outlier if its
  \code{\link{ks.test}} test statistic is more than \code{theta}
  standard deviations larger than the test statistics.}
 }

\value{An integer vector of indices (see \code{\link[base:which]{which}}).}

\details{
The function is very simple, have a look at its code.
}

\author{Wolfgang Huber <whuber@embl.de> }

\seealso{\code{\link[stats:ks.test]{ks.test}}.}

