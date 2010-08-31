\name{aqm.prepaffy}
\alias{aqm.prepaffy}

\title{ Preparation of Affymetrix experiments for RLE and NUSE. }

\description{
  \code{aqm.prepaffy} performs data preprocessing on
  \code{\link[affy:AffyBatch-class]{AffyBatch}} and formats the output
  for \code{aqm.rle} and \code{aqm.nuse} usage.
}

\usage{ aqm.prepaffy(expressionset) }

\arguments{
  \item{expressionset}{ is an object of
    class \code{\link[affy:AffyBatch-class]{AffyBatch}} }
}

\value{ A preprocessed affy object of class
  \code{\link{aqmobj.prepaffy}}. }

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqm.rle}}, \code{\link{aqm.nuse}} }

\examples{
    library(ALLMLL)
    data(MLL.A)
    MLLprep = aqm.prepdata(MLL.A, do.logtransform = TRUE)
    MLLaffyprep = aqm.prepaffy(MLL.A)
    nuse = aqm.nuse(MLL.A, MLLprep, MLLaffyprep)
    class(nuse)
    aqm.plot(nuse)
}