## Spatial plot representation
spatialplot = function(expressionset, dataprep, channel, label, imageMat = NULL)
  {
    noaxis = function (...) 
      {     
        tick="no"
      }                        
    
    colourRamp = colorRampPalette(rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256)))
    if(dataprep$classori != "AffyBatch")
      {
        r = featureData(expressionset)$X
        c = featureData(expressionset)$Y
      }
    if(dataprep$classori == "AffyBatch")
      {
        maxc = ncol(expressionset)
        maxr = nrow(expressionset)

        r = rep(seq_len(maxr), maxc)
        c = rep(seq_len(maxc), each = maxr)
      }
    if(dataprep$classori == "BeadLevelList")
      {        
        maxc = ncol(imageMat[[1]])
        maxr = nrow(imageMat[[1]])
        
        r = rep(seq_len(maxr), maxc)
        c = rep(seq_len(maxc), each = maxr)
      }
    df = if(dataprep$classori == "BeadLevelList") data.frame("Array" = rep(seq_len(dataprep$numArrays), each=maxr*maxc), "ch" = unlist(imageMat), "row" = r, "column" = c) else data.frame("Array" = as.factor(col(dataprep$dat)),  "ch" = unlist(lapply(seq_len(dataprep$numArrays), function(i) rank(channel[,i]))), "row" = r,  "column" = c)
    levelplot(ch ~ column*row | Array, data=df, axis = noaxis, asp = "iso", col.regions = colourRamp, as.table=TRUE, strip = function(..., bg) strip.default(..., bg ="#cce6ff"), xlab = label, ylab = "") #
  }

## Scores computation
scoresspat = function(expressionset, dataprep, ch)
  {
    if(dataprep$classori == "BeadLevelList")
      stop("Cannot apply scoresspat to BeadLevelList objects. Please use scoresspatbll instead.")
    if(dataprep$classori != "AffyBatch")
      {
        maxr = max(featureData(expressionset)$X)
        maxc = max(featureData(expressionset)$Y)
      }
    if(dataprep$classori == "AffyBatch")
      {
        maxr = nrow(expressionset)
        maxc = ncol(expressionset)
      }
    
    mdat = lapply(seq_len(dataprep$numArrays), function(x) matrix(ch[,x],ncol=maxc,nrow=maxr,byrow=T))
    loc = sapply(mdat, function(x) {
      apg = abs(fft(x)) ## absolute values of the periodogram
      lowFreq = apg[1:4, 1:4]
      lowFreq[1,1] = 0                  # drop the constant component
      highFreq = c(apg[-(1:4), ], apg[1:4, -(1:4)])
      return(sum(lowFreq)/sum(highFreq))
    })
    return(loc)
  }


scoresspatbll = function(expressionset, dataprep, log)
{
  if(dataprep$classori == "BeadLevelList")
  {        
    imageMatg = lapply(seq_len(dataprep$numArrays), function(i) imageplotBLL(expressionset, array = i, whatToPlot="G", log = log))
    if(expressionset@arrayInfo$channels == "two")
      imageMatr = lapply(seq_len(dataprep$numArrays), function(i) imageplotBLL(expressionset, array = i, whatToPlot="Rb", log = log))
    
    locr = sapply(imageMatr, function(x) {
      x[x == "NaN"] = 0
      apg = abs(fft(x)) ## absolute values of the periodogram
      lowFreq = apg[1:4, 1:4]
      lowFreq[1,1] = 0  # drop the constant component
      highFreq = c(apg[-(1:4), ], apg[1:4, -(1:4)])
      return(sum(lowFreq)/sum(highFreq))
    })
     locg = sapply(imageMatg, function(x) {
      x[x == "NaN"] = 0
      apg = abs(fft(x)) ## absolute values of the periodogram
      lowFreq = apg[1:4, 1:4]
      lowFreq[1,1] = 0  # drop the constant component
      highFreq = c(apg[-(1:4), ], apg[1:4, -(1:4)])
      return(sum(lowFreq)/sum(highFreq))
    })
    loc = list("Red" = locr, "Green" = locg)
    return(list(locr, locg))
  }
}

## set Generics  
setGeneric("spatialbg",
           function(expressionset, dataprep)
           standardGeneric("spatialbg"))
setGeneric("spatial",
           function(expressionset, dataprep)
           standardGeneric("spatial"))

## NCS
##Background rank representation
setMethod("spatialbg",signature(expressionset = "NChannelSet"), function(expressionset, dataprep)
          {
                sprb = spatialplot(expressionset, dataprep,dataprep$rcb,"Rank(red background intensity)")
                spgb = spatialplot(expressionset, dataprep,dataprep$gcb,"Rank(green background intensity)")
                return(list(sprb, spgb))
          })
            
## NCS
##Foreground rank representation
setMethod("spatial",signature(expressionset = "NChannelSet"), function(expressionset, dataprep)
          {
            spr = spatialplot(expressionset, dataprep,dataprep$rc,"Rank(red intensity)")
            spg = spatialplot(expressionset, dataprep,dataprep$gc,"Rank(green intensity)")
            splr = spatialplot(expressionset, dataprep,dataprep$dat,"Rank(log(ratio))")
            return(list(spr, spg, splr))
          })
              
## ES-AB
## Foreground rank representation
setMethod("spatial",signature(expressionset = "aqmOneCol"), function(expressionset, dataprep)
          {
            spi = spatialplot(expressionset, dataprep,dataprep$dat,"Rank(intensity)")
            return(spi)
          })
            
##BeadLeveList plot modified from beadarry package
imageplotBLL = function(BLData, array = 1, nrow = 100, ncol = 100,
                     whatToPlot ="G", log=TRUE, zlim=NULL, method="illumina",
                     n = 3, trim=0.05, ...){

  whatToPlot = match.arg(whatToPlot, choices=c("G", "Gb", "R", "Rb", "wtsG", "wtsR", "residG", "residR", "M", "residM", "A", "beta"))
  if((whatToPlot=="R" | whatToPlot=="residR" | whatToPlot=="M" | whatToPlot=="residM" | whatToPlot=="A" | whatToPlot=="beta") & BLData@arrayInfo$channels!="two")
    stop(paste("Need two-channel data to plot", whatToPlot, "values"))
                                          
  data = getArrayData(BLData, what=whatToPlot, array=array, log=log, method=method, n=n, trim=trim) 
  ind = is.na(data) | is.infinite(data)
  if(sum(ind)>0) {
    cat(paste("Warning:", sum(ind), "NA, NaN or Inf values, which will be ignored.\nCheck your data or try setting log=\"FALSE\"\n"))
    data = data[!ind]
  }

  xs = floor(BLData[[array]]$GrnX[!ind])
  ys = floor(BLData[[array]]$GrnY[!ind])
  rm(ind)

  xgrid = floor(seq(0, max(xs), by = max(xs)/ncol))
  ygrid = floor(seq(0, max(ys), by = max(ys)/nrow))

  imageMatrix = matrix(ncol = ncol, nrow = nrow)
  zr = NULL
  for(i in seq_len(ncol)){
    idx = which((xs > xgrid[i]) & (xs < xgrid[i+1]))
    if(length(idx)>0) {
      fground = data[idx]
      yvalues = ys[idx]

      out = .C("BLImagePlot", length(fground), as.double(fground), as.double(yvalues), as.integer(ygrid),
        result = double(length = nrow), as.integer(nrow), PACKAGE = "beadarray")
      zr[1] = min(zr[1], out$result[!is.na(out$result)], na.rm=TRUE)
      zr[2] = max(zr[2], out$result[!is.na(out$result)], na.rm=TRUE)
      if(!is.null(zlim)) {
        out$result[!is.na(out$result)] = pmax(zlim[1], out$result[!is.na(out$result)], na.rm=TRUE)
        out$result[!is.na(out$result)] = pmin(zlim[2], out$result[!is.na(out$result)], na.rm=TRUE)
      }
      imageMatrix[,i] = rev(out$result)
    }
  }
  return(imageMatrix)
}

## BLL
setMethod("spatialbg",signature(expressionset = "BeadLevelList"), function(expressionset, dataprep)
          {
            imageMatgb = lapply(seq_len(dataprep$numArrays), function(i) imageplotBLL(expressionset, array = i, whatToPlot="Gb", log = dataprep$logtransformed))
            if(expressionset@arrayInfo$channels == "two")
              imageMatrb = lapply(seq_len(dataprep$numArrays), function(i) imageplotBLL(expressionset, array = i, whatToPlot="Rb", log = dataprep$logtransformed))
            
            spgb = spatialplot(expressionset, dataprep, channel, label, imageMatgb)
            sprb = NULL
            if(expressionset@arrayInfo$channels == "two")
              sprb = spatialplot(expressionset, dataprep, channel, label, imageMatrb)
            
            return(list(sprb, spgb))
          })
            
setMethod("spatial",signature(expressionset = "BeadLevelList"), function(expressionset, dataprep)
          {
            imageMatg = lapply(seq_len(dataprep$numArrays), function(i) imageplotBLL(expressionset, array = i, whatToPlot="G", log = dataprep$logtransformed))
            if(expressionset@arrayInfo$channels == "two")
              imageMatr = lapply(seq_len(dataprep$numArrays), function(i) imageplotBLL(expressionset, array = i, whatToPlot="R", log = dataprep$logtransformed))
            
            spg = spatialplot(expressionset, dataprep, channel, label, imageMatg)
            spr = NULL
            if(expressionset@arrayInfo$channels == "two")
              spr = spatialplot(expressionset, dataprep, channel, label, imageMatr)
            
            return(list(spr, spg))
          })

## BACKGROUND
aqm.spatialbg = function(expressionset, dataprep)
{

  bg = spatialbg(expressionset, dataprep)
  
  title = "Spatial distribution of local background intensities"
  type = "Individual array quality"
  
  legend = sprintf("<b>Figure FIG:</b> False color representations of the arrays' spatial distributions of feature local background estimates. The color scale is shown in the panel on the right, and it is proportional to the ranks of the probe background intensities.")

  out = list("plot" = bg, "type" = type, "title" = title, "legend" = legend)
  class(out) = "aqmobj.spatialbg"
  return(out)
}

## FOREGROUND
aqm.spatial = function(expressionset, dataprep)
{
  fore = spatial(expressionset, dataprep)
  title = "Spatial distribution of feature intensities"
  type = "Individual array quality"

  legend = sprintf("<b>Figure FIG:</b> False color representations of the arrays' spatial distributions of feature intensities. The color scale is shown in the panel on the right, and it is proportional to the ranks of the probe intensities. This has indeed the potential to detect patterns that are small in amplitude, but systematic within array. These may not be consequential for the downstream data analysis, but if properly interpreted, could e.g. still be useful for technology and experimental protocol optimisation as it helps in identifying patterns that may be caused by, for example, spatial gradients in the hybridization chamber, air bubbles, spotting or plating problems.")          

  if(!is(expressionset, "BeadLevelList"))
    {
      loc = scoresspat(expressionset, dataprep, dataprep$dat)
  
      if(is(expressionset, "NChannelSet"))
        {
          locr = scoresspat(expressionset, dataprep, dataprep$rc)
          locg = scoresspat(expressionset, dataprep, dataprep$gc)
          loc = list("LogRatio"=loc, "Red"=locr, "Green"=locg)
        }

      if(is(loc,"numeric"))
        {
          locstat = boxplot.stats(loc)
          locout = sapply(seq_len(length(locstat$out)), function(x) which(loc == locstat$out[x]))
        }
 
      if(is(loc,"list"))
        {
          locstat = lapply(seq_len(length(loc)), function(i) boxplot.stats(loc[[i]]))
          locout = lapply(seq_len(length(loc)), function(j) sapply(seq_len(length(locstat[[j]]$out)), function(x) which(loc[[j]] == locstat[[j]]$out[x])))
          names(locout) = c("LR","Red","Green")
        }
    }

  if(is(expressionset, "BeadLevelList"))
    loc = scoresspatbll(expressionset, dataprep, dataprep$logtransformed)

  out = list("plot" = fore, "type" = type, "title" = title, "legend" = legend, "scores" = loc, "outliers" = locout)
  class(out) = "aqmobj.spatial"
  return(out)
}
