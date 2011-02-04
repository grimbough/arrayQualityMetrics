library("arrayQualityMetrics")
options(error=recover)


## Download a set of Human-WG6 V2 Expression Arrays
## from Mike Smith's webpage
## http://www.compbio.group.cam.ac.uk/Resources/baloc
arrayNames = paste(
  rep(c("4343238066", "4343238080"), each=12),
  rep(LETTERS[1:6], each=2, times=2),
  rep(c("1", "2"), times=12),
  sep = "_")

for(f in paste(arrayNames, "tar.gz", sep=".")) {
  if(!file.exists(f)) {
    download.file(paste("http://www.compbio.group.cam.ac.uk/Resources/baloc/Data", f, sep="/"),
                  destfile=f)
  }
}

##
## Uncompress
##
library("BeadDataPackR")
for(f in arrayNames) {
  f1 = paste(f, "tar.gz", sep=".")
  f2 = paste(f, "bab", sep=".")
  f3 = paste(f, "Grn.tif", sep=".")
  f4 = paste(f, "txt", sep=".")
  if(!file.exists(f4)) {
    system(paste("tar -xzvf", f1))
    decompressBeadData(input=f2)
    file.remove(c(f2, f3))
  }
}

##
## read them in
##
library("beadarray")
if(!exists("BLData")) {
  if(!file.exists("BLData.rda")) {
    BLData = readIllumina(illuminaAnnotation = "Humanv2")
    save(BLData, file="BLData.rda") 
  } else {
    load("BLData.rda")
  }
}


## work around an apparent bug in beadarray:
BLData@sectionData$Metrics = NULL

##
## summarise
##
if(!exists("Bexpr")) {
  if(!file.exists("Bexpr.rda")) {
    
    myMean = function(x) mean(x, na.rm = TRUE)
    mySd = function(x) sd(x, na.rm = TRUE)
    green = new("illuminaChannel", logGreenChannelTransform, illuminaOutlierMethod, myMean, mySd, "G")
    Bexpr = summarize(BLData, list(green), useSampleFac = FALSE)
    save(Bexpr, file="Bexpr.rda")
    
  } else {
    load("Bexpr.rda")

  }
}

arrayQualityMetrics(Bexpr, force=TRUE)

