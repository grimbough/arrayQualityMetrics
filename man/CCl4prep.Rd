\name{CCl4prep}
\docType{CCl4prep}
\alias{CCl4prep}
\title{aqmobj.prepdata object from an example experiment}
\description{
 \code{CCl4prep} is an object of class \code{aqmobj.prepdata} which contains 
  the CCl4norm data prepared by the function \code{aqm.prepdata}.}
\usage{data(CCl4prep)}
\seealso{\code{\link[arrayQualityMetrics:class-aqmobj.prepdata]{aqmobj.prepdata}}}
\examples{
    data(CCl4norm)
    CCl4prep = aqm.prepdata(CCl4norm, do.logtransform = FALSE)
}
\keyword{datasets}
