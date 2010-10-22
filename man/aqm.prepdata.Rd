\name{prepdata}
\alias{prepdata}
\alias{prepaffy}

\title{Compute useful summary statistics from a data object.}

\usage{
prepdata(expressionset, intgroup, do.logtransform)
prepaffy(expressionset)
}

\arguments{
  \item{expressionset}{An object of class
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} for one colour
  non Affymetrix data,
  \code{\link[affy:AffyBatch-class]{AffyBatch}} for Affymetrix data,
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}} for two colour arrays, or
  \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}} for Illumina
  bead arrays.}
  \item{do.logtransform}{Logical of length 1, indicating whether to log
    transform the data.}
}

\description{\code{prepdata} computes summary statistics that are useful
  for all platforms; \code{prepaffy} Affymetrix-specific ones.}

\value{A list with various derived values.}
\author{Audrey Kauffmann.}

\examples{
    ## Load an example of a NChannelSet
    library("CCl4")
    data("CCl4")

    ## Normalization of CCl4 using vsn
    library("vsn")
    CCl4norm = justvsn(CCl4, subsample=2000)

    ## Add a column in the phenoData to annotate samples 
    cond = paste(CCl4norm$RIN.Cy3, CCl4norm$RIN.Cy5, sep="/")
    poor   = grep(cond, pattern="2.5")
    medium = grep(cond, pattern="^5/|/5")
    good   = grep(cond, pattern="9.7")
    cov = rep(0, length = nrow(pData(CCl4norm)))
    cov[good]   = "Good"
    cov[medium] = "Medium"
    cov[poor]   = "Poor"
    phenoData(CCl4norm)$RNAintegrity = cov

    ## Add X and Y columns in the featureData to allow spatial representations
    featureData(CCl4norm)$X = featureData(CCl4norm)$Row
    featureData(CCl4norm)$Y = featureData(CCl4norm)$Column

    ## Add a 'hasTarget' column in the featureData, to enable the
    ## aqm.probesmap plot
    featureData(CCl4norm)$hasTarget = (regexpr("^NM",
                                       featureData(CCl4norm)$Name)> 0)

    ## Call prepdata
    CCl4prep = prepdata(CCl4norm, do.logtransform = FALSE, intgroup = "RNAintegrity")

    ## Draw MA plots
    ma = aqm.maplot(CCl4prep)
    class(ma)
    aqm.plot(ma)

    ## Draw heatmap making use of the RNAintegrity
    ## column of the phenoData
    hm = aqm.heatmap(CCl4prep)
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
    pm = aqm.probesmap(CCl4prep)
    class(pm)
    aqm.plot(pm)
}
