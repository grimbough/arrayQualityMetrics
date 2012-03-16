#------------------------------------------------------------
# Tests on RGList and MAList objects
#------------------------------------------------------------
library("arrayQualityMetrics")
library("limma")
options(error=recover)

what = c("RGList", "MAList")

if("RGList" %in% what){
  library("CCl4")

  data("CCl4_RGList")
  data("CCl4")
  stopifnot( is.data.frame(CCl4_RGList$targets),
             identical( paste(rownames(pData(CCl4)), ".gpr", sep=""), CCl4_RGList$targets$FileName ))
  CCl4_RGList$targets = pData(CCl4)

  arrayQualityMetrics(CCl4_RGList, intgroup="RIN.Cy3", do.logtransform=TRUE, force=TRUE)
}

if("MAList" %in% what){
  ## From http://www-huber.embl.de/users/whuber/bioc-list/120315
  load("MA.avg_TEST1.RData")
  arrayQualityMetrics(MA.avg_TEST1, intgroup="Cy3", force=TRUE)
}

