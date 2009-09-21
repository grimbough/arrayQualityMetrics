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
  \code{\link[affy:AffyBatch-class]{AffyBatch}}, a
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}}, or a
  \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}} into a
  \code{\link{aqmobj.prepdata}} object which can be used as an input of 
  the aqm functions.}

\usage{
aqm.prepdata(expressionset, do.logtransform = TRUE)
}

\arguments{
  \item{expressionset}{ An object of class
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} for one colour
  non Affymetrix data, \code{\link[affy:AffyBatch-class]{AffyBatch}} for
  Affymetrix data, \code{\link[Biobase:class.NChannelSet]{NChannelSet}}
  for two colour arrays, or
  \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}} for Illumina
  bead arrays.
 }
  \item{do.logtransform}{ TRUE or FALSE whether or not you want to log
  transform the data.}
 }

\value{ An object of class  \code{\link{aqmobj.prepdata}}. }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\seealso{\code{\link{aqmobj.prepdata}}, 
\code{\link{aqm.boxplot}}, 
\code{\link{aqm.density}}, 
\code{\link{aqm.heatmap}}, 
\code{\link{aqm.maplot}}, 
\code{\link{aqm.meansd}}, 
\code{\link{aqm.probesmap}}, 
\code{\link{aqm.spatial}}, 
\code{\link{aqm.spatialbg}}}

\examples{
    ## Load an example of a NChannelSet
    library(CCl4)
    data(CCl4)

    ## Normalization of CCl4 using vsn
    library(vsn)
    CCl4norm = justvsn(CCl4, subsample=2000)

    ## Add a column in the phenoData to annotate samples 
    cond = paste(pData(CCl4norm)$RIN.Cy3,pData(CCl4norm)$RIN.Cy5,sep="/")
    poor = grep(cond,pattern="2.5")
    medium = grep(cond,pattern="^5/|/5")
    good = grep(cond,pattern="9.7")
    cov = rep(0, length = nrow(pData(CCl4norm)))
    cov[good] = "Good"
    cov[medium] = "Medium"
    cov[poor] = "Poor"
    phenoData(CCl4norm)$RNAintegrity = cov

    ## Add X and Y columns in the featureData to allow spatial representations
    featureData(CCl4norm)$X = featureData(CCl4norm)$Row
    featureData(CCl4norm)$Y = featureData(CCl4norm)$Column

    ## Add a hasTarget column in the featureData to call aqm.probesmap
    featureData(CCl4norm)$hasTarget = (regexpr("^NM",
                                       featureData(CCl4norm)$Name)> 0)


    ## Prepare the data for aqm.xxx calls
    CCl4prep = aqm.prepdata(CCl4norm, do.logtransform = FALSE)

    ## Draw MA plots
    ma = aqm.maplot(dataprep = CCl4prep)
    class(ma)
    aqm.plot(ma)

    ## Draw heatmap making use of the RNAintegrity
    ## column of the phenoData
    hm = aqm.heatmap(expressionset = CCl4norm,
                     dataprep = CCl4prep,
                     intgroup = "RNAintegrity")
    class(hm)
    aqm.plot(hm)

    ## Draw probes mapping density curves making use of the hasTarget
    ## column of the featureData
    sp = aqm.spatial(expressionset = CCl4norm,
                     dataprep = CCl4prep,
                     scale = "Rank")
    class(sp)
    aqm.plot(sp)

    ## Draw probes mapping density curves making use of the hasTarget
    ## column of the featureData
    pm = aqm.probesmap(expressionset = CCl4norm, dataprep = CCl4prep)
    class(pm)
    aqm.plot(pm)


}