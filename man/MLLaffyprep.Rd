\name{MLLaffyprep}
\docType{MLLaffyprep}
\alias{MLLaffyprep}
\title{aqmobj.prepaffy object from an example experiment}
\description{
 \code{MLLaffyprep} is an object of class \code{aqmobj.prepaffy} which contains 
  the MLL.A data prepared by the function \code{aqm.prepaffy}.}
\usage{data(MLLaffyprep)}
\seealso{\code{\link[arrayQualityMetrics:class-aqmobj.prepaffy]{aqmobj.prepaffy}}}
\examples{
    library(ALLMLL)
    data(MLL.A)
    MLLaffyprep = aqm.prepaffy(MLL.A, sampleNames(MLL.A))
}
\keyword{datasets}
