\name{aqm.qcstats}
\alias{aqm.qcstats}

\title{ Performs QCstats plot on AffyBatch. }

\description{
  \code{aqm.qcstats} performs QCstats on objects of class
  \code{\link[affy:AffyBatch-class]{AffyBatch}} and formats the output for
  \code{aqm.plot} usage.
}

\usage{ aqm.qcstats(expressionset, ...) }

\arguments{
  \item{expressionset}{ is an object of
    class \code{\link[affy:AffyBatch-class]{AffyBatch}} }
  \item{...}{ Any arguments to \code{\link[simpleaffy:qc]{qc}} }
 }

\value{ An object of class  \code{\link{aqmobj.qcs}}. }

\author{ Audrey Kauffmann <kauffmann@bergonie.org> }

\seealso{ \code{\link{aqmobj.qcs}} }

\examples{
    library(ALLMLL)
    data(MLL.A)
    qm = aqm.qcstats(MLL.A)
    class(qm)
    aqm.plot(qm)
}
