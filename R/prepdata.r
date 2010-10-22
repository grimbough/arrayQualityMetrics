setClassUnion("aqmOneCol", c("ExpressionSet", "AffyBatch"))

logtransform = function(x) {
  if (is.null(x)) {
    NULL
  } else { 
    ifelse(x>0, suppressWarnings(log2(x)), rep(NA_real_, length(x)))
  }
}

## Platform specific data-preparation
setGeneric("platformspecific",
           function(expressionset,
                    do.logtransform = TRUE)
           standardGeneric("platformspecific"))

## NChannelSset
setMethod("platformspecific",
          signature(expressionset = "NChannelSet"),
function(expressionset, do.logtransform){
  rc = assayData(expressionset)$R
  gc = assayData(expressionset)$G
  
  rcb = if(exists("Rb", envir = assayData(expressionset))) assayData(expressionset)$Rb else NULL
  gcb = if(exists("Gb", envir = assayData(expressionset))) assayData(expressionset)$Gb else NULL

  if(do.logtransform)
    {
      rc  = logtransform(rc)
      gc  = logtransform(gc)
      rcb = logtransform(rcb)
      gcb = logtransform(gcb)
    } 

  M = dat = rc - gc
  A = 0.5*(rc+gc)
  
  colnames(dat) = seq(along = sampleNames(expressionset))
  numArrays = ncol(rc)
  
  if("dyeswap" %in% names(phenoData(expressionset)@data))
    {
      lev = levels(expressionset@phenoData$dyeswap)
      if(length(lev) != 2)
        stop("The dyeswap slot of the phenoData must be binary.\n")
      reverseddye = names(expressionset@phenoData$dyeswap[expressionset@phenoData$dyeswap == min(lev)])
      dat[,reverseddye] = - dat[,reverseddye]
    }
  
  outM = as.dist(dist2(dat))
  
  list(
     "rc" = rc, "gc" = gc, "rcb" = rcb, "gcb" = gcb,
     "M" = M, "A" = A, "dat" = dat,
     "outM" = outM, "numArrays" = numArrays,
     "nchannels" = 2)
})

## ExpressionSet and AffyBatch
setMethod("platformspecific",
          signature(expressionset = "aqmOneCol"),
function(expressionset, do.logtransform){
  dat = exprs(expressionset)    
  if(do.logtransform)
    dat = logtransform(dat)
  
  medArray = rowMedians(dat, na.rm=TRUE)
  M =  dat - medArray
  A =  (dat + medArray)/2
  
  numArrays = ncol(dat)
  outM = as.dist(dist2(dat))
  
  list(
    "rc" = NULL, "gc" = NULL,  "rcb" = NULL, "gcb" = NULL,
    "M" = M, "A" = A, "dat" = dat,
    "outM" = outM, "numArrays" = numArrays,
    "nchannels" = 1)
})
          
## BeadLevelList
setMethod("platformspecific",
          signature(expressionset = "BeadLevelList"),
function(expressionset, do.logtransform){
            
  numArrays = as.numeric(dim(expressionset)[1])

  switch(expressionset@arrayInfo$channels,
   "single" =  {
     summaryES = createBeadSummaryData(expressionset, imagesPerArray = 1, log = do.logtransform)
     rc = gc = NULL
     nch = 1
   },
   "two" = {
     summaryESRG = createBeadSummaryData(expressionset, what = "RG", imagesPerArray = 1, log = do.logtransform)
     rc = assayData(summaryESRG)$R                         
     gc = assayData(summaryESRG)$G                
     summaryES = createBeadSummaryData(expressionset ,what = "M", imagesPerArray = 1, log = do.logtransform)
     nch = 2
   },
   stop(sprintf("Invalid expressionset@arrayInfo$channels = '%s'", expressionset@arrayInfo$channels))    
  ) ## switch
  
  dat = exprs(summaryES)
  outM = as.dist(dist2(dat))
 
  switch(expressionset@arrayInfo$channels,
   "single" = {
      medArray = rowMedians(dat, na.rm=TRUE)
      M =  dat - medArray
      A =  (dat + medArray)/2
   },
   "two" = {
      M = rc - gc
      A = 0.5*(rc +gc)
   })
  
  list(
    "rc" = rc, "gc" = gc, "rcb" = NULL, "gcb" = NULL,
    "M" = M, "A" = A, "dat" = dat,
    "outM" = outM, "numArrays" = numArrays,
    "nchannels" = nch)
})


prepdata = function(expressionset, intgroup, do.logtransform) {

  cls = class(expressionset)
  if (any(cls %in% c("RGList", "MAList", "marrayRaw", "marrayNorm"))) {
    expressionset = try(as(expressionset, "NChannelSet"))
    if(inherits(expressionset,'try-error'))
      stop(sprintf("Argument 'expressionset' is of class '%s', and its automatic conversion into 'NChannelSet' failed. Please try to convert it manually.\n", paste(cls, collapse=", ")))
  }

  x = platformspecific(expressionset, do.logtransform)

  x$intgroup = intgroup
  x$do.logtransform = do.logtransform
  x$classori = class(expressionset)[1]
  x$expressionset = expressionset
  
  return(x)
}
