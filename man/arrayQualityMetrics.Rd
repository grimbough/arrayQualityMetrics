\name{arrayQualityMetrics-package}
\alias{arrayQualityMetrics}
\docType{package}
\title{
Quality metrics on ExpressionSets
}
\description{
  This package performs quality metrics on ExpressionSet containing 
  microarray data from any platforms, one or two channels. 
  The results are designated to allow the user to rapidly assess the
  quality of a set of arrays.
}
\usage{
arrayQualityMetrics(expressionset, outfile, do.logtransform, split.plots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expressionset}{ The \code{expressionset} must be an object of
    class ExpressionSet }
  \item{outfile}{ The name of the result file is \code{outfile} }
  \item{do.logtransform}{ default is FALSE. If TRUE, the data will be
    log transformed }
  \item{split.plots}{  default is FALSE. If the number of studied array
    is more than 50 it is adviced to define a number of experiments to
    represent on the density plots }
}
\details{
\tabular{ll}{
Package: \tab arrayQualityMetrics\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-06-29\cr
License: \tab What license is it under?\cr
}
ATTENTION An overview of how to use the package, including the most important 
 functions 
}
\author{
Audrey Kauffmann, Wolfgang Huber
Maintainer: <audrey@ebi.ac.uk>
}
\references{
  ATTENTION ref
}
\keyword{ package }
\seealso{
See Also \code{\link[<pkg>:<pkg>-package]{<pkg>}} 
}
\examples{
}
\keyword{ hplot }
\keyword{ graphs }% __ONLY ONE__ keyword per line
