\name{aqm.qcstats}

\alias{aqm.qcstats}
\alias{aqm.pmmm}
\alias{aqm.rnadeg}
\alias{aqm.rle}
\alias{aqm.nuse}

\title{Compute Affymetrix specific plots.}

\description{
  \code{aqm.qcstats} performs QCstats on objects of class
  \code{\link[affy:AffyBatch-class]{AffyBatch}} and formats the output for
  \code{aqm.plot} usage.
}

\usage{
aqm.qcstats(expressionset)
aqm.pmmm(expressionset)
aqm.rnadeg(expressionset)
aqm.rle(x, affyproc)
aqm.nuse(x, affyproc)
}

\arguments{
  \item{expressionset}{An \code{AffyBatch} object.}
  \item{x}{The object resulting from the call to \code{\link{prepdata}(expressionset)}.}
  \item{affyproc}{The object resulting from the call to \code{\link{prepaffy}(expressionset)}.}
}

\value{
An object of class \code{aqmTrellis} or \code{aqmPlotfun}.
}

\author{Audrey Kauffmann}

\seealso{\code{\link{aqm.boxplot}}}

\examples{
    library("ALLMLL")
    data(MLL.A)
    qm = aqm.qcstats(MLL.A)
    class(qm)
    aqm.plot(qm)
}
