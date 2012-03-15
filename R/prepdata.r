prepdata = function(expressionset, intgroup, do.logtransform) {
  conversions = c(`RGList` = "NChannelSet")
  for(i in seq(along=conversions)) {
    if (is(expressionset, names(conversions)[i])) {
      expressionset = try(as(expressionset, conversions[i]))
      if(is(expressionset, "try-error")) {
        stop(sprintf("The argument 'expressionset' is of class '%s', and its automatic conversion into '%s' failed. Please try to convert it manually, or contact the creator of that object.\n", names(conversions)[i], conversions[i]))
      } else {
        break
      }
    }
  }

  x = platformspecific(expressionset, do.logtransform)

  if (!is.null(intgroup)) {
    if (!is.character(intgroup))
      stop("'intgroup' should be a 'character'.")
    if(!all(intgroup %in% colnames(x$pData)))
      stop("all elements of 'intgroup' should match column names of 'pData(expressionset)'.")
  }

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

##--------------------------------------------------
## Platform specific data-preparation
##--------------------------------------------------
setGeneric("platformspecific",
           function(expressionset,
                    do.logtransform)
           standardGeneric("platformspecific"))

##--------------------------------------------------
## NChannelSset (either one or two colors)
##--------------------------------------------------
setMethod("platformspecific",
          signature(expressionset = "NChannelSet"),
function(expressionset, do.logtransform) {

  adNames = ls(assayData(expressionset))
  nIsOne = ("exprs" %in% adNames)
  nIsTwo = all(c("R","G") %in% adNames)
  if(!xor(nIsOne, nIsTwo))
      error("'assayData(expressionset)' must contain either 'exprs', or 'R' and 'G', but not both.")

  if(nIsOne)
      return(oneColor(expressionset, do.logtransform=do.logtransform, M=assayData(expressionset)$exprs))

  if(!nIsTwo)
      error("'assayData(expressionset)' must contain either 'exprs', or 'R' and 'G'.")

  adNotUsed = setdiff(adNames, c("R", "G", "Rb", "Gb"))
  if(length(adNotUsed)>0)
      warning(sprintf("Elements%s %s of 'assayData(expressionset)' will be ignored.",
                      if(length(adNotUsed)==1) "" else "s",
                      paste(adNotUsed, collapse=", ")))

  R  = assayData(expressionset)$R   ## red channel foreground
  G  = assayData(expressionset)$G   ## green channel foreground
  Rb = assayData(expressionset)$Rb  ## red channel background
  Gb = assayData(expressionset)$Gb  ## green channel background
  sx = featureData(expressionset)$X ## spatial x-coordinate
  sy = featureData(expressionset)$Y ## spatial y-coordinate

  if(do.logtransform) {
    R  = logtransform(R)
    G  = logtransform(G)
    Rb = logtransform(Rb)
    Gb = logtransform(Gb)
  }

  M = R-G
  A = 0.5*(R+G)

  M = applyDyeSwap(M, phenoData(expressionset))

  list(R = R, G = G, Rb = Rb, Gb = Gb,
       sx = sx, sy = sy,
       M = M, A = A,
       nchannels = 2, pData = cleanPhenoData(expressionset), fData = fData(expressionset))
})


##--------------------------------------------------
## MAList
##--------------------------------------------------
setMethod("platformspecific",
          signature(expressionset = "MAList"),
function(expressionset, do.logtransform)
{

  if(!identical(FALSE, do.logtransform))
    warning("Ignoring argument 'do.logtransform', since there is no such choice for 'MAList' objects.")

  M = expressionset$M
  A = expressionset$A
  R = A+M/2
  G = A-M/2

  pd = cleanUpPhenoData(expressionset$targets)
  M = applyDyeSwap(M, pd)

  fd = expressionset$genes

  list(R = R, G = G, Rb = NULL, Gb = NULL,
       sx = NULL, sy = NULL,
       M = M, A = A,
       nchannels = 2, pData = pd, fData = fd)
})

applyDyeSwap = function(M, pd) {
  if("dyeswap" %in% colnames(pd)) {
    if(!is.factor(pd$dyeswap))
        stop("'pd$dyeswap' must be a factor.")
    lev = levels(pd$dyeswap)
    if(length(lev) != 2)
        stop("The factor 'pd$dyeswap' must have exactly two levels")

    rev = as.integer(pd$dyeswap)==1
    M[, rev] = -M[, rev]
  }
  return(M)
}

##----------------------------------------------------------
## Common parts of dealing with ExpressionSet and AffyBatch
##----------------------------------------------------------
oneColor = function(expressionset, do.logtransform, M)
{
  if(missing(M))
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
function(expressionset, do.logtransform) {
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
function(expressionset, do.logtransform)
    oneColor(expressionset, do.logtransform))

##------------------------------------------------------------
## extract and clean up phenoData
##------------------------------------------------------------
cleanPhenoData = function(x, ...) {

  pd = pData(x)

  scd = protocolData(x)[["ScanDate"]]
  if(!is.null(scd))
    if(length(scd)==nrow(pd))
      pd = cbind(pd, ScanDate=scd)

  cleanUpPhenoData(pd, ...)
}

cleanUpPhenoData = function(pd, maxcol = 10) {

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
