##--------------------------------------------------------------------------------
## Spatial distributions. This function calls 'spatialplot' for possibly multiple
## elements of the list x, including "M", "R", "G"
##--------------------------------------------------------------------------------
aqm.spatial = function(x, scale="rank", channels = c("M", "R", "G"))
{
  if(!( (length(scale)==1) && is.character(scale) && (scale %in% c("direct", "rank")) ))
    stop("'scale' must be 'direct' or 'rank'\n")

  channels = intersect(channels, names(x))
  
  if(is.numeric(x$sx) && is.numeric(x$sy))
    {
      lapply(channels, spatialplot, x=x, scale=scale)
    } else {
      list(new("aqmReportModule",
          plot = NULL,
          section = "Individual array quality",
          title = paste("Spatial intensity distribution", if(length(channels)>1) "s" else "", sep=""),
          legend = "Spatial distributions could not be plotted. Please provide columns <tt>sx</tt> and <tt>sy</tt> with the features' spatial coordinates in the <tt>featureData</tt> slot of the data object."))
    }
}


##--------------------------------------------------------------------------------
## Spatial distribution plot for one channel 
##--------------------------------------------------------------------------------
spatialplot = function(whichChannel, x, scale) 
{
    
  colourRamp = colorRampPalette(rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256)))

  dat = x[[whichChannel]]
  stopifnot(length(x$sx)==nrow(dat),
            length(x$sy)==nrow(dat))  ## this should always be true, given the definition of prepdata
  
  if(scale=="rank")
    dat = apply(dat, 2, rank)
  
  df = data.frame(
    "Array"  = as.factor(col(dat)),
    "ch"     = as.vector(dat),
    "row"    = x$sy,
    "column" = x$sx)
  
  spat = levelplot(ch ~ column*row | Array,
            data = df,
            axis = function (...) list(tick="no"),
            asp = "iso",
            col.regions = colourRamp,
            as.table = TRUE,
            strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
            xlab = whichChannel,
            ylab = "",
            panel = "panel.levelplot.raster",
            colorkey = (scale!="rank"))
         
  legend = sprintf("The figure <!-- FIG --> shows false colour representations of the arrays' spatial distributions of feature intensities. Normally, when the features are distributed randomly on the arrays, one expects to see a uniform distribution; sets of control features with particularly high or low intensities may stand out. The colour scale is proportional to %sthe probe intensities, and it is shown in the panel on the right.", switch(scale, rank = "the ranks of ", direct = ""))

  if(scale=="rank") legend = paste(legend, "Note that the rank scale has the potential to amplify patterns that are small in amplitude but systematic within an array. It is possible to switch off the rank scaling by modifying the argument 'scale' in the call of the 'aqm.spatial' function.") 
  
  ## Outlier detection
  loc = scoresspat(x, whichChannel)
  locstat = boxplot.stats(loc)
  locout  = which(loc %in% locstat$out)
  
  new("aqmReportModule",
      plot = spat,
      section = "Individual array quality",
      title = paste("Spatial distribution of", whichChannel),
      legend = legend,
      outliers = locout,
      shape = list("h"=6,"w"=6))
}

  
## Scores computation
scoresspat = function(x, wh)
{

  loc = numeric(x$numArrays)
  maxx = max(x$sx, na.rm=TRUE)
  maxy = max(x$sy, na.rm=TRUE)
  
  for(a in seq_len(x$numArrays))
    {
      mat = matrix(NA_real_, nrow=maxy, ncol=maxx)
      mat[ cbind(x$sy, x$sx) ] = x[[wh]][, a]
      apg = abs(fft(mat))       ## absolute values of the periodogram
      lowFreq      = apg[1:4, 1:4]
      lowFreq[1,1] = 0          ## drop the constant component
      highFreq     = c(apg[-(1:4), ], apg[1:4, -(1:4)])
      loc[a] = sum(lowFreq)/sum(highFreq)
    }
  return(loc)
}



