\name{aqm.nuse}
\alias{aqm.nuse}

\title{ Performs NUSE plot on aqmobj.prepaffy objects. }

\description{
  \code{aqm.nuse} performs NUSE boxplots and outlier detection from it and
  formats the output for \code{aqm.plot} usage.  
}

\usage{ aqm.nuse(affyproc, ...) }

\arguments{
  \item{affyproc}{ An object of class \code{\link{aqmobj.prepaffy}} }
  \item{...}{ Any arguments to \code{boxplot} }
}

\value{ An object of class  \code{\link{aqmobj.nuse}} }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{ \code{\link{aqm.prepaffy}}, \code{\link{aqmobj.prepaffy}},
  \code{\link{aqmobj.nuse}} }

\examples{
    library(ALLMLL)
    data(MLL.A)
    pa = aqm.prepaffy(MLL.A, sampleNames(MLL.A))
    nuse = aqm.nuse(pa)
}