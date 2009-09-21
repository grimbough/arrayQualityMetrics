\name{aqm.pmmm}
\alias{aqm.pmmm}

\title{ Performs perfect match versus mismatch density plots. }

\description{
  \code{aqm.pmmm} performs PM MM density curves on objects of class
  \code{\link[affy:AffyBatch-class]{AffyBatch}} and formats the output for
  \code{aqm.plot} usage.
}

\usage{ aqm.pmmm(expressionset, ...) }

\arguments{
  \item{expressionset}{ is an object of
    class \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{...}{ Any arguments to \code{density} }
}

\value{ An object of class  \code{\link{aqmobj.pmmm}}. }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqmobj.pmmm}} }

\examples{
    library(ALLMLL)
    data(MLL.A)
    pm = aqm.pmmm(MLL.A)
    class(pm)
    aqm.plot(pm)
}
