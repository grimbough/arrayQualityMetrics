\name{aqm.plot}
\alias{aqm.plot}
\alias{aqm.plot-methods}
\alias{aqm.plot,aqmTrellis-method}
\alias{aqm.plot,aqmobj.box-method}
\alias{aqm.plot,aqmobj.dens-method}
\alias{aqm.plot,aqmobj.msd-method}
\alias{aqm.plot,aqmobj.nuse-method}
\alias{aqm.plot,aqmobj.pmmm-method}
\alias{aqm.plot,aqmobj.qcs-method}
\alias{aqm.plot,aqmobj.rle-method}
\alias{aqm.plot,aqmobj.rnadeg-method}

\title{ Performs plots from aqm objects. }

\description{ \code{aqm.plot} performs plots. }

\usage{
\S4method{aqm.plot}{aqmTrellis}(obj)

aqm.plot(obj)

\S4method{aqm.plot}{aqmobj.box}(obj)
\S4method{aqm.plot}{aqmobj.dens}(obj)
\S4method{aqm.plot}{aqmobj.msd}(obj)
\S4method{aqm.plot}{aqmobj.nuse}(obj)
\S4method{aqm.plot}{aqmobj.pmmm}(obj)
\S4method{aqm.plot}{aqmobj.qcs}(obj)
\S4method{aqm.plot}{aqmobj.rle}(obj)
}

\arguments{
  \item{obj}{ an object of class \code{\link[arrayQualityMetrics:class.aqmobj]{aqmobj}}}
 }

\value{ A plot. }

\author{
  Audrey Kauffmann. 
  Maintainer: <audrey@ebi.ac.uk>  
}

\examples{
  \dontrun{
    library(ALLMLL)
    data(MLL.A)
    dataprep = aqm.prepdata(MLL.A, TRUE)
    ma = aqm.maplot(dataprep)
    aqm.plot(ma)
    }
}