prepdata = function(expressionset, intgroup, do.logtransform)
{
  cls = class(expressionset)
  if (is(expressionset, "RGList") || is(expressionset, "MAList"))
    {
      expressionset = try(as(expressionset, "NChannelSet"))
      if(is(expressionset, "try-error"))
        stop(sprintf("Argument 'expressionset' is of class '%s', and its automatic conversion into 'NChannelSet' failed. Please try to convert it manually.\n", paste(cls, collapse=", ")))
    }

  x = platformspecific(expressionset, do.logtransform)  # see below

  x = append(x, list(
    numArrays       = ncol(x$M),
    intgroup        = intgroup,
    do.logtransform = do.logtransform))

  x = append(x, intgroupColors(x))
  
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
       nchannels = 2, pData = cleanPhenoData(expressionset), fData = fData(expressionset))
})

##----------------------------------------------------------
## Common parts of dealing with ExpressionSet and AffyBatch
##----------------------------------------------------------
oneColor = function(expressionset, do.logtransform)
{
  M = exprs(expressionset)    
  if(do.logtransform)
    M = logtransform(M)
  
  list(M = M, A = M, 
       nchannels = 1, pData = cleanPhenoData(expressionset), fData = fData(expressionset))
}

##----------------------------------------------------------
## ExpressionSet
##----------------------------------------------------------
setMethod("platformspecific",
          signature(expressionset = "ExpressionSet"),
function(expressionset, do.logtransform)
{
  rv = oneColor(expressionset, do.logtransform)
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
  rv = oneColor(expressionset, do.logtransform)
  maxc = ncol(expressionset)
  maxr = nrow(expressionset)
  rv$sx = rep(seq_len(maxc), each = maxr) ## spatial x-coordinate
  rv$sy = rep(seq_len(maxr), maxc)        ## spatial y-coordinate
  return(rv)
})

##----------------------------------------------------------
## beadLevelData 
##----------------------------------------------------------
setMethod("platformspecific",
          signature(expressionset = "beadLevelData"),
function(expressionset, do.logtransform)
{
  stop("\n\narrayQualityMetrics does not support bead-level objects of class 'beadLevelData'. Please refer to the vignette 'Analysis of bead-level data using beadarray' of the beadarray package for the processing and quality assessment of such objects. After summarisation of the bead-level data to an object of class 'ExpressionSetIllumina', you can use arrayQualityMetrics on that.\n")
})

##----------------------------------------------------------
## ExpressionSetIllumina
##----------------------------------------------------------
setMethod("platformspecific",
          signature(expressionset = "ExpressionSetIllumina"),
          oneColor)

##------------------------------------------------------------
## extract and clean up phenoData
##------------------------------------------------------------
cleanPhenoData = function(x, maxcol = 10)
{

  pd = pData(x)

  scd = protocolData(x)[["ScanDate"]]
  if(!is.null(scd))
    if(length(scd)==nrow(pd))
      pd = cbind(pd, ScanDate=scd)
  
  if(ncol(pd) > maxcol)
    {
      ## Remove columns whose contents are all the same, or all different (except for the first time).
      ## After that, if there are still more than maxcol columns left, only use the first ones.
      remove = rep(FALSE, ncol(pd))
      hadUniqueAlready = FALSE
      for(i in seq_len(ncol(pd)))
        {
          n = nlevels(factor(pd[[i]]))
          if(n==1)
            {
              remove[i] = TRUE
            } else {
              if(n==nrow(pd))
                {
                  if(hadUniqueAlready)
                    remove[i] = TRUE else hadUniqueAlready = TRUE
                }
            }
        } ## for i
      
      wh = which(!remove)
      stopifnot(length(wh)>0)
      if(length(wh) > maxcol)
        wh = wh[seq_len(maxcol)]
      
      pd = pd[, wh, drop=FALSE]
    } ## if ncol(x)
  
  ## Convert everything which is not a factor into a character string, and remove "@"
  for(i in seq_len(ncol(pd)))
    {
      if(is.factor(pd[[i]]))
        {
          levels(pd[[i]]) = gsub("@", " ", levels(pd[[i]]), fixed=TRUE)
        } else {
          pd[[i]] = gsub("@", " ", as.character(pd[[i]]), fixed=TRUE)
        }
    }
  
  return(pd)
}
