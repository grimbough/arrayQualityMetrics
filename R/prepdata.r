logtransform = function(x) {
  if (is.null(x)) {
    NULL
  } else {
    pos = (x>0)
    x[!pos] = NA_real_
    x[pos] = log2(x[pos])
    x
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
  
  rc  = assayData(expressionset)$R   ## red channel foreground
  gc  = assayData(expressionset)$G   ## green channel foreground
  rcb = assayData(expressionset)$Rb  ## red channel background
  gcb = assayData(expressionset)$Gb  ## green channel background
  sx  = featureData(expressionset)$X ## spatial x-coordinate
  sy  = featureData(expressionset)$Y ## spatial y-coordinate

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
  
  if("dyeswap" %in% colnames(phenoData(expressionset)))
    {
      if(!is.factor(expressionset$dyeswap))
        stop("'expressionset$dyeswap' must be a factor.")
      lev = levels(expressionset$dyeswap)
      if(length(lev) != 2)
        stop("The factor 'expressionset$dyeswap' must have exactly two levels")

      reverse = (as.integer(expressionset$dyeswap)==1)
      dat[,rev] = - dat[,rev]
    }
  
  list(
     "rc" = rc, "gc" = gc, "rcb" = rcb, "gcb" = gcb, sx = sx, sy = sy,
     "M" = M, "A" = A, "dat" = dat,
     "nchannels" = 2, pData = pData(expressionset), fData = fData(expressionset))
})

## ExpressionSet and AffyBatch
setMethod("platformspecific",
          signature(expressionset = "oneColourArray"),
function(expressionset, do.logtransform){
  
  dat = exprs(expressionset)    
  if(do.logtransform)
    dat = logtransform(dat)
  
  medArray = rowMedians(dat, na.rm=TRUE)
  M =  dat - medArray
  A =  (dat + medArray)/2
  
  list(
    "M" = M, "A" = A, "dat" = dat,
    "nchannels" = 1, pData = pData(expressionset), fData = fData(expressionset))
})
          
## BeadLevelList
setMethod("platformspecific",
          signature(expressionset = "BeadLevelList"),
function(expressionset, do.logtransform){
            
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
     summaryES = createBeadSummaryData(expressionset , what = "M", imagesPerArray = 1, log = do.logtransform)
     nch = 2
   },
   stop(sprintf("Invalid 'expressionset@arrayInfo$channels': %s", expressionset@arrayInfo$channels))    
  ) ## switch
  
  dat = exprs(summaryES)
 
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
    "rc" = rc, "gc" = gc, 
    "M" = M, "A" = A, "dat" = dat,
    "nchannels" = nch, pData = pData(expressionset), fData = fData(expressionset)) ## TODO - this needs to be fixed.
})


prepdata = function(expressionset, intgroup, do.logtransform, usesvg) {

  cls = class(expressionset)
  if (any(cls %in% c("RGList", "MAList", "marrayRaw", "marrayNorm"))) {
    expressionset = try(as(expressionset, "NChannelSet"))
    if(inherits(expressionset,'try-error'))
      stop(sprintf("Argument 'expressionset' is of class '%s', and its automatic conversion into 'NChannelSet' failed. Please try to convert it manually.\n", paste(cls, collapse=", ")))
  }

  x = platformspecific(expressionset, do.logtransform)

  x = append(x, list(
    outM            = as.dist(dist2(x$dat)),
    numArrays       = ncol(x$dat),
    intgroup        = intgroup,
    do.logtransform = do.logtransform,
    usesvg          = usesvg,
    classori        = class(expressionset)[1],  ## TODO is this needed anywhere? If not, drop
    expressionset   = expressionset))
  
  return(x)
}
