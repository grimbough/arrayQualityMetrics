## Methods to prepare the data
prepaffy = function(expressionset)
{
  pp1 = preprocess(expressionset)
  dataPLM = fitPLM(pp1, background = FALSE, normalize = FALSE)
  list("dataPLM" = dataPLM)
}

## phrase book
phrase = list(
  fig = function(x) sprintf("The figure <!-- FIG --> shows the <i>%s</i> plot.", x),
  preproc = function(x) sprintf("%s values are computed from the preprocessed data (after background correction and quantile normalisation).", x),
  affyPLM = "This diagnostic plot is based on tools provided in the <i>affyPLM</i> package."
  )

## RNA digestion
aqm.rnadeg = function(expressionset, x)
  {
 
    cl = intgroupColours(x)
    
    rnaDeg = function() {
      plotAffyRNAdeg(AffyRNAdeg(expressionset, log.it = TRUE),
                     lwd = 1, col = cl$arrayColours)
    }
    
    legend = paste(phrase$fig("RNA digestion"),
                   phrase$preproc("Shown"),
       "In this plot each array is represented by a single line; move the mouse over the lines to see their corresponding sample names. The plot can be used to identify array(s) that have a slope very different from the others. This could indicate that the RNA used for that array has been handled differently from what was done for the other arrays.",
                   phrase$affyPLM)

    new("aqmReportModule",
        plot    = rnaDeg,
        section = "Affymetrix specific plots",
        title   = "RNA digestion plot",
        legend  = legend,
        shape   = list("h" = 5.5, "w" =7),
        svg     =  if(x$usesvg)
        list(getfun = aqm.getMatplotSeries, numObjects = x$numArrays,
             strokewidth = c(1,3), strokeopacity = c(0.4, 1)) else list())
   }


## RLE
aqm.rle = function(x, affyproc)
  {   
    if (affyproc$dataPLM@model.description$R.model$which.parameter.types[3] == 1){
      medianchip = apply(coefs(affyproc$dataPLM), 1, median)
      dat = sweep(coefs(affyproc$dataPLM), 1, medianchip, FUN='-')
    } else {
      stop("It doesn't appear that a model with sample effects was used.")
    }

    rle = aqm.boxplot(list(dat=dat, intgroup=x$intgroup, nchannels=1, numArrays=x$numArrays))

    rle@legend = paste(phrase$fig("Relative Log Expression (RLE)"),
      phrase$preproc("RLE"),
      "Arrays whose boxes have larger spread or are centered away from M = 0 are potentially problematic.",
      phrase$affyPLM)
    
    rle@title = "Relative Log Expression plot"
    rle@section = "Affymetrix specific plots"
    rle@outliers = integer(0)

    ## TODO: Fix the outlier computation
    ## rle2 = try(Mbox(affyproc$dataPLM, plot = FALSE))
    ## rlemed = rle2$stats[3,]
    ## rleout = which(abs(rlemed) > 0.1)

    return(rle)
  }

## NUSE
aqm.nuse = function(x, affyproc)
{ 
  ## bwplot for PLMset
  ## TODO: Use 'colSums' - it's faster
  compute.nuse <- function(which){
    1/sqrt(apply(affyproc$dataPLM@weights[[1]][which,,drop=FALSE], 2, sum))
  }
  
  model <- affyproc$dataPLM@model.description$modelsettings$model
  if ((affyproc$dataPLM@model.description$R.model$which.parameter.types[3] == 1) &&
      (affyproc$dataPLM@model.description$R.model$which.parameter.types[1] == 0) ){
    grp.rma.se1.median  = apply(se(affyproc$dataPLM), 1,median,na.rm=TRUE)
    grp.rma.rel.se1.mtx = sweep(se(affyproc$dataPLM),1,grp.rma.se1.median,FUN='/')
  } else {
    which <- indexProbesProcessed(affyproc$dataPLM)
    ses   <- matrix(0, length(which), 4)
    if (affyproc$dataPLM@model.description$R.model$response.variable == 1){
      for (i in seq(along=which))
        ses[i,] <- compute.nuse(which[[i]])
    } else {
      stop("Sorry, I don't know how to compute NUSE values for this PLMset object")
    }
    grp.rma.se1.median <- apply(ses, 1,median)
    grp.rma.rel.se1.mtx <- sweep(ses,1,grp.rma.se1.median,FUN='/')
  }

  nuse = aqm.boxplot(list(dat=grp.rma.rel.se1.mtx, intgroup=x$intgroup, nchannels=1, numArrays=x$numArrays))

  nuse@legend = paste(phrase$fig("Normalized Unscaled Standard Error (NUSE)"),
    phrase$preproc("NUSE"),
    "Potential low quality arrays are those whose box is substantially elevated or more spread out relative to the other arrays. NUSE values are comparable within one data set, but usually not across different data sets.",
    phrase$affyPLM)

  nuse@title = "Normalized Unscaled Standard Error (NUSE) plot"
  nuse@section = "Affymetrix specific plots"
  nuse@outliers = integer(0)

  return(nuse)
}


## QCStats
aqm.qcstats = function(expressionset) {                

  qcStats = function() {
    plot.qc.stats(qc(expressionset))
  }
  
  legend =  paste(phrase$fig("the Affymetrix recommended diagnostic"),
    "Please see the vignette of the package <i>simpleaffy</i> for a full explanation of the elements shown in this plot. Any metrics that is shown in red is outside the manufacturer's specified boundaries and suggests a potential problem; metrics shown in blue are considered acceptable.")
  
  shape = list("h" = 4 + ncol(exprs(expressionset)) * 0.1 + 1/ncol(exprs(expressionset)),  "w" = 6)
  
  new("aqmReportModule",
      plot = qcStats,
      section = "Affymetrix specific plots",
      title = "Diagnostic plot recommended by Affymetrix",
      legend = legend,
      shape = shape)
}

## PM / MM
aqm.pmmm = function(expressionset)
{  
  PM = density(as.matrix(log2(pm(expressionset))))
  MM = density(as.matrix(log2(mm(expressionset))))

  PMMM = function(){
    plot(MM, col = "grey", xlab = "log(Intensity)", main="")
    lines(PM, col = "blue")
    legend("topright", c("PM","MM"), lty=1, lwd=2, col= c("blue","grey"), bty = "n")
  }
  
  legend = "Figure <!-- FIG --> shows the density distributions of the log<sub>2</sub> intensities grouped by the matching type of the probes. The blue line shows a density estimate (smoothed histogram) from intensities of perfect match probes (PM), the grey line, one from the mismatch probes (MM). We expect that MM probes have poorer hybridization than PM probes, and thus that the PM curve be to the right of the MM curve."

  new("aqmReportModule",
      plot = PMMM,
      section = "Affymetrix specific plots",
      title = "Perfect matches and mismatches",
      legend = legend,
      shape = list("h" = 6, "w" = 6))
}



