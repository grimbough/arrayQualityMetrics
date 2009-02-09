\name{aqm.prepdata}
\alias{aqm.prepdata}
\alias{aqm.prepdata,ExpressionSet-method}
\alias{aqm.prepdata,AffyBatch-method}
\alias{aqm.prepdata,NChannelSet-method}
\alias{aqm.prepdata,BeadLevelList-method}
\alias{aqm.prepdata,aqmOneCol-method}

\title{ Generate an object aqmobj.prepdata to be called by the aqm functions. }

\description{ \code{aqm.prepdata} formats an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}, an
  \code{\link[affy:class.AffyBatch]{AffyBatch}}, a
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}}, or a
  \code{\link[beadarray:class.BeadLevelList]{BeadLevelList}} into a
  \code{\link{aqmobj.prepdata}} object which can be used as an input of
  the aqm functions.}

\usage{
\S4method{aqm.prepdata}{ExpressionSet}(expressionset, do.logtransform)

aqm.prepdata(expressionset, do.logtransform = TRUE)

\S4method{aqm.prepdata}{AffyBatch}(expressionset, do.logtransform)
\S4method{aqm.prepdata}{NChannelSet}(expressionset, do.logtransform)
\S4method{aqm.prepdata}{BeadLevelList}(expressionset, do.logtransform)
\S4method{aqm.prepdata}{aqmOneCol}(expressionset, do.logtransform)
}

\arguments{
  \item{expressionset}{ An object of class
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} for one colour
  non Affymetrix data, \code{\link[affy:class.AffyBatch]{AffyBatch}} for
  Affymetrix data, \code{\link[Biobase:class.NChannelSet]{NChannelSet}}
  for two colour arrays, or
  \code{\link[beadarray:class.BeadLevelList]{BeadLevelList}} for Illumina
  bead arrays.
 }
  \item{do.logtransform}{ TRUE or FALSE whether or not you want to log
  transform the data.}
 }

\value{ An object of class  \code{\link{aqmobj.prepdata}}. }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{\code{\link{aqmobj.prepdata}}, \code{\link{aqm.boxplot}}, \code{\link{aqm.density}}, \code{\link{aqm.heatmap}}, \code{\link{aqm.maplot}}, \code{\link{aqm.meansd}}, \code{\link{aqm.probesmap}}, \code{\link{aqm.spatial}}, \code{\link{aqm.spatialbg}}}

\examples{
    data(CCl4norm)
    ## Prepare the data for aqm.xxx calls
    CCl4prep = aqm.prepdata(CCl4norm, do.logtransform = FALSE)
}