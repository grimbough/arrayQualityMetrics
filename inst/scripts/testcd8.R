##
## Test a dataset from Tim Rayner (see his email of 10 Nov 2010)
## cd8 is an AffyBatch object created from HuGene ST 1.0
##

library("arrayQualityMetrics")
options(error=recover)

if(!file.exists("cd8.Rdata")) 
  download.file("http://dl.dropbox.com/u/1225281/cd8.RData",
                  destfile="cd8.Rdata")

if(!exists("cd8"))
  load("cd8.Rdata")

arrayQualityMetrics(cd8, force=TRUE, do.logtransform=TRUE)

