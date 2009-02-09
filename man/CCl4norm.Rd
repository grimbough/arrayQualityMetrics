\name{CCl4norm}
\docType{CCl4norm}
\alias{CCl4norm}
\title{NChannelSet object from an example experiment}
\description{
 \code{CCl4norm} is an object of class \code{NChannelSet} which contains 
  the CCl4 data normalised and with additionnal annotation in the featureData and phenoData.}
\usage{data(CCl4norm)}
\seealso{\code{\link[Biobase:class-NChannelSet]{NChannelSet}}}
\examples{
    ## Load CCl4 data
    library(CCl4)
    library(vsn)
    data(CCl4)

    ## Add a column in the phenoData to indicate the integrity of the RNA
    p = pData(CCl4)
    cond = paste(p$RIN.Cy3,p$RIN.Cy5,sep="/")
    poor = grep(cond,pattern="2.5")
    medium = grep(cond,pattern="^5/|/5")
    good = grep(cond,pattern="9.7")
    cov = rep(0, length = nrow(p))
    cov[good] = "Good"
    cov[medium] = "Medium"
    cov[poor] = "Poor"
    phenoData(CCl4)$RNAintegrity = cov

    ## Call columns Row and Column from the featureData X and Y
    fData(CCl4)$X = fData(CCl4)$Row
    fData(CCl4)$Y = fData(CCl4)$Column

    ## Add a column hasTarget containing TRUE and FALSE if the probe
    ## name starts with 'NM' or not
    fData(CCl4)$hasTarget = (regexpr("^NM", fData(CCl4)$Name) > 0)

    ## Normalise with VSN
    CCl4norm = justvsn(CCl4,subsample=2200)
}
\keyword{datasets}
