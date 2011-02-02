prepdata = function(expressionset, intgroup, do.logtransform)
{
  cls = class(expressionset)
  if (any(cls %in% c("RGList", "MAList", "marrayRaw", "marrayNorm")))
    {
      expressionset = try(as(expressionset, "NChannelSet"))
      if(inherits(expressionset, "try-error"))
        stop(sprintf("Argument 'expressionset' is of class '%s', and its automatic conversion into 'NChannelSet' failed. Please try to convert it manually.\n", paste(cls, collapse=", ")))
    }
  
  x = platformspecific(expressionset, do.logtransform)  # see below

  x = append(x, list(
    numArrays       = ncol(x$M),
    intgroup        = intgroup,
    do.logtransform = do.logtransform))
  
  return(x)
}


##--------------------------------------------------
## Little helper function
##--------------------------------------------------
logtransform = function(x)
{
  if (is.null(x)) {
    NULL
  } else {
    pos = (x>0)
    x[!pos] = NA_real_
    x[pos] = log2(x[pos])
    x
  }
}

##--------------------------------------------------
## Platform specific data-preparation
##--------------------------------------------------
setGeneric("platformspecific",
           function(expressionset,
                    do.logtransform = TRUE)
           standardGeneric("platformspecific"))

##--------------------------------------------------
## NChannelSset (two colors)
##--------------------------------------------------
setMethod("platformspecific",
          signature(expressionset = "NChannelSet"),
function(expressionset, do.logtransform)
{
  R  = assayData(expressionset)$R   ## red channel foreground
  G  = assayData(expressionset)$G   ## green channel foreground
  Rb = assayData(expressionset)$Rb  ## red channel background
  Gb = assayData(expressionset)$Gb  ## green channel background
  sx = featureData(expressionset)$X ## spatial x-coordinate
  sy = featureData(expressionset)$Y ## spatial y-coordinate

  if(do.logtransform)
    {
      R  = logtransform(R)
      G  = logtransform(G)
      Rb = logtransform(Rb)
      Gb = logtransform(Gb)
    } 

  M = R-G
  A = 0.5*(R+G)
  
  if("dyeswap" %in% colnames(phenoData(expressionset)))
    {
      if(!is.factor(expressionset$dyeswap))
        stop("'expressionset$dyeswap' must be a factor.")
      lev = levels(expressionset$dyeswap)
      if(length(lev) != 2)
        stop("The factor 'expressionset$dyeswap' must have exactly two levels")

      rev = as.integer(expressionset$dyeswap)==1
      M[,rev] = -M[,rev]
    }
  
  list(R = R, G = G, Rb = Rb, Gb = Gb,
       sx = sx, sy = sy,
       M = M, A = A, 
       nchannels = 2, pData = pData(expressionset), fData = fData(expressionset))
})

##----------------------------------------------------------
## Common parts of dealing with ExpressionSet and AffyBatch
##----------------------------------------------------------
oneColourStuff = function(expressionset, do.logtransform)
{
  M = exprs(expressionset)    
  if(do.logtransform)
    M = logtransform(M)
  
  list(M = M, A = M, 
       nchannels = 1, pData = pData(expressionset), fData = fData(expressionset))
}

##----------------------------------------------------------
## ExpressionSet
##----------------------------------------------------------
setMethod("platformspecific",
          signature(expressionset = "ExpressionSet"),
function(expressionset, do.logtransform)
{
  rv = oneColourStuff(expressionset, do.logtransform)
  rv$sx = featureData(expressionset)$X ## spatial x-coordinate
  rv$sy = featureData(expressionset)$Y ## spatial y-coordinate
  return(rv)
})
          
##----------------------------------------------------------
## AffyBatch
##----------------------------------------------------------
setMethod("platformspecific",
          signature(expressionset = "AffyBatch"),
function(expressionset, do.logtransform)
{
  rv = oneColourStuff(expressionset, do.logtransform)

  maxc = ncol(expressionset)
  maxr = nrow(expressionset)
  rv$sx = rep(seq_len(maxc), each = maxr) ## spatial x-coordinate
  rv$sy = rep(seq_len(maxr), maxc)        ## spatial y-coordinate
  return(rv)
})

##----------------------------------------------------------
## BeadLevelList
## TODO - this needs to be fixed; also need to extract x and y coordinates (sx, sy)
##----------------------------------------------------------
setMethod("platformspecific",
          signature(expressionset = "BeadLevelList"),
function(expressionset, do.logtransform){
            
  switch(expressionset@arrayInfo$channels,
   "single" =  {
     summaryES = createBeadSummaryData(expressionset, imagesPerArray = 1, log = do.logtransform)
     R = G = NULL
     M = A = exprs(summaryES) 
     nch = 1
   },
   "two" = {
     summaryESRG = createBeadSummaryData(expressionset, what = "RG", imagesPerArray = 1, log = do.logtransform)
     R = assayData(summaryESRG)$R                         
     G = assayData(summaryESRG)$G                
     summaryES = createBeadSummaryData(expressionset , what = "M", imagesPerArray = 1, log = do.logtransform)
     M = exprs(summaryES)   
     A = 0.5*(R +G)         ## TODO us there also a createBeadSummaryData function for this?
     nch = 2
   },
   stop(sprintf("Invalid 'expressionset@arrayInfo$channels': %s", expressionset@arrayInfo$channels))    
  ) ## switch
  
  list(R = R, G = G, M = M, A = A, 
    nchannels = nch, pData = pData(expressionset), fData = fData(expressionset)) 
})


