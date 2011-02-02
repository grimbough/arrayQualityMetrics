##--------------------------------------------------
## Phrases
##--------------------------------------------------
figurePhrase = function(x)
  sprintf("The figure <!-- FIG --> shows the <i>%s</i> plot.", x)

##--------------------------------------------------
## Prepare the data
##--------------------------------------------------
prepaffy = function(expressionset, x)
{
  x$dataPLM = fitPLM(expressionset)
  x$pm      = pm(expressionset)
  x$mm      = mm(expressionset)
  x$mmOK    = !any(is.na(x$mm))
  return(x)
}

##--------------------------------------------------
## RNA digestion
##--------------------------------------------------
aqm.rnadeg = function(expressionset, x)
  {
 
    cl = intgroupColours(x)
    
    rnaDeg = function() {
      plotAffyRNAdeg(AffyRNAdeg(expressionset, log.it = TRUE),
                     lwd = 1, col = cl$arrayColours)
    }
    
    legend = paste(figurePhrase("RNA digestion"),
      "The shown values are computed from the preprocessed data (after background correction and quantile normalisation). Each array is represented by a single line; move the mouse over the lines to see their corresponding sample names. The plot can be used to identify array(s) that have a slope very different from the others. This could indicate that the RNA used for that array has been handled differently from what was done for the other arrays.")

    new("aqmReportModule",
        plot    = rnaDeg,
        section = "Affymetrix specific plots",
        title   = "RNA digestion plot",
        legend  = legend,
        shape   = list("h" = 5.5, "w" =7),
        svg     =  new("svgParameters",
                       name = "rnadeg",
                       numPlotObjects = x$numArrays) )

   }

##--------------------------------------------------
## RLE
##--------------------------------------------------
aqm.rle = function(x)
{
  outlierMethod = "KS"
  values = RLE(x$dataPLM, type="values")
  rv = aqm.boxplot(list(M=values, intgroup=x$intgroup, nchannels=1, numArrays=x$numArrays), outlierMethod=outlierMethod)

  rv@title = "Relative Log Expression (RLE)"
  rv@section = "Affymetrix specific plots"
  
  rv@legend = paste(figurePhrase(rv@title),
    "Arrays whose boxes are centered away from 0 and/or are more spread out are potentially problematic.",
    outlierPhrase(outlierMethod, length(rv@outliers)))

  return(rv)
}

##--------------------------------------------------
## NUSE
##--------------------------------------------------
aqm.nuse = function(x)
{
  outlierMethod = "median"
  values = NUSE(x$dataPLM, type="values")
  rv = aqm.boxplot(list(M=values, intgroup=x$intgroup, nchannels=1, numArrays=x$numArrays), outlierMethod=outlierMethod)

  rv@title = "Normalized Unscaled Standard Error (NUSE)"
  rv@section = "Affymetrix specific plots"

  rv@legend = paste(figurePhrase(rv@title),
    "For each array, the boxes should be centered around 1. An array were the values elevated relative to the other arrays is typically of lower quality.",
    outlierPhrase(outlierMethod, length(rv@outliers)))

  return(rv)
}

##--------------------------------------------------
## simpleaffy 
## If fail, silently return NULL
##--------------------------------------------------
## aqm.qcstats = function(expressionset) {                

##   qcObj = try(qc(expressionset))

##   if(inherits(qcObj, "try-error"))
##     return(NULL)
  
##   qcStats = function() 
##     plot.qc.stats(qcObj)
    
##   legend =  paste(figurePhrase("diagnostics suggested by Affymetrix"),
##     "Please see the vignette of the package <i>simpleaffy</i> for a full explanation of the elements shown in this plot. Any metrics that is shown in red is outside the manufacturer's specified boundaries and suggests a potential problem; metrics shown in blue are considered acceptable.")
  
##   shape = list("h" = 4 + ncol(exprs(expressionset)) * 0.1 + 1/ncol(exprs(expressionset)),  "w" = 6)
  
##   new("aqmReportModule",
##       plot = qcStats,
##       section = "Affymetrix specific plots",
##       title = "Diagnostic plot recommended by Affymetrix",
##       legend = legend,
##       shape = shape)
## }

##--------------------------------------------------
## PM / MM
## If fail, silently return NULL
##--------------------------------------------------
aqm.pmmm = function(x)
{
  if(!x$mmOK)
    return(NULL)

  PM = density(as.matrix(log2(x$pm)))
  MM = density(as.matrix(log2(x$mm)))
  
  PMMM = function(){
    plot(MM, col = "grey", xlab = "log(Intensity)", main="")
    lines(PM, col = "blue")
    legend("topright", c("PM", "MM"), lty=1, lwd=2, col=c("blue","grey"), bty="n")
  }
  
  legend = "Figure <!-- FIG --> shows the density distributions of the log<sub>2</sub> intensities grouped by the matching type of the probes. The blue line shows a density estimate (smoothed histogram) from intensities of perfect match probes (PM), the grey line, one from the mismatch probes (MM). We expect that MM probes have poorer hybridization than PM probes, and thus that the PM curve be to the right of the MM curve."
  
  new("aqmReportModule",
      plot = PMMM,
      section = "Affymetrix specific plots",
      title = "Perfect matches and mismatches",
      legend = legend,
      shape = list("h" = 6, "w" = 6))
  
}



