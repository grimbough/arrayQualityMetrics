setClassUnion("aqmOneCol", c("ExpressionSet", "AffyBatch"))

##Functions
logtransform = function(x)
  ifelse(x>0, suppressWarnings(log2(x)), rep(NA_real_, length(x)))

##Methods to prepare the data
setGeneric("aqm.prepdata",
           function(expressionset,
                    do.logtransform = TRUE)
           standardGeneric("aqm.prepdata"))

##NCS
setMethod("aqm.prepdata",signature(expressionset = "NChannelSet"), function(expressionset, do.logtransform){
  if(do.logtransform)
    {
      rc = logtransform(assayData(expressionset)$R)
      gc = logtransform(assayData(expressionset)$G)
    } else {
      rc = assayData(expressionset)$R
      gc = assayData(expressionset)$G
    }
  if("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))
    {
      rcb = logtransform(assayData(expressionset)$Rb)
      gcb = logtransform(assayData(expressionset)$Gb)
    } else rcb = gcb = NULL

  M = dat = rc - gc
  A = 0.5*(rc+gc)
  
  sN = seq_len(length(sampleNames(expressionset)))
  
  numArrays = ncol(rc)
  colnames(dat) = sN
  if("dyeswap" %in% names(phenoData(expressionset)@data))
    {
      lev = levels(expressionset@phenoData$dyeswap)
      if(length(lev) != 2)
        stop("The dyeswap slot of the phenoData must be binary.\n")
      reverseddye = names(expressionset@phenoData$dyeswap[expressionset@phenoData$dyeswap == min(lev)])
      dat[,reverseddye] = - dat[,reverseddye]
    }
  outM = as.dist(dist2(dat))
  
  
  object = list("rc" = rc, "gc" = gc, "rcb" = rcb, "gcb" = gcb, "M" = M, "A" = A, "dat" = dat, "outM" = outM, "sN" = sN, "numArrays" = numArrays, "nchannels" = 2, "logtransformed" = do.logtransform, "classori" = class(expressionset)[1])
  class(object) = "aqmobj.prepdata"
  return(object)
})

##ES & AB
setMethod("aqm.prepdata",signature(expressionset = "aqmOneCol"), function(expressionset, do.logtransform){
  if(do.logtransform)
    {
      dat = logtransform(exprs(expressionset))
    } else dat = exprs(expressionset)    
  
  medArray = rowMedians(dat, na.rm=TRUE)
  M =  dat - medArray
  A =  (dat + medArray)/2
  
  sN = seq_len(length(sampleNames(expressionset)))
 
  numArrays = ncol(dat)
  outM = as.dist(dist2(dat))
  
  object = list("rc" = NULL, "gc" = NULL,  "rcb" = NULL, "gcb" = NULL,"M" = M, "A" = A, "dat" = dat, "outM" = outM, "sN" = sN, "numArrays" = numArrays, "nchannels" = 1, "logtransformed" = do.logtransform, "classori" = class(expressionset)[1])
  class(object) = "aqmobj.prepdata"
  return(object)
})
          
##BLL
setMethod("aqm.prepdata",signature(expressionset = "BeadLevelList"), function(expressionset, do.logtransform){
            
 sN = seq_len(length(arrayNames(expressionset)))
 
  numArrays = as.numeric(dim(expressionset)[1])
            
  if(expressionset@arrayInfo$channels == "single")
    {
      summaryES = createBeadSummaryData(expressionset, imagesPerArray = 1, log = do.logtransform)
      rc = gc = NULL
      nch = 1
    }
  if(expressionset@arrayInfo$channels == "two")
    {
      summaryESRG = createBeadSummaryData(expressionset, what = "RG", imagesPerArray = 1, log = do.logtransform)
      rc = assayData(summaryESRG)$R                         
      gc = assayData(summaryESRG)$G                
      summaryES = createBeadSummaryData(expressionset ,what = "M", imagesPerArray = 1, log = do.logtransform)
      nch = 2
    }            
  
  dat = exprs(summaryES)
  outM = as.dist(dist2(dat))
 
  if(expressionset@arrayInfo$channels == "single")
    {
      medArray = rowMedians(dat, na.rm=TRUE)
      M =  dat - medArray
      A =  (dat + medArray)/2
    }
  if(expressionset@arrayInfo$channels == "two")
    {
      M = rc - gc
      A = 0.5*(rc +gc)
    }
  
  object = list("rc" = rc, "gc" = gc, "rcb" = NULL, "gcb" = NULL, "M" = M, "A" = A, "dat" = dat, "outM" = outM, "sN" = sN, "numArrays" = numArrays, "nchannels" = nch, "logtransformed" = do.logtransform, "classori" = class(expressionset)[1])
  class(object) = "aqmobj.prepdata"
  return(object)
})

