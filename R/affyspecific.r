##Methods to prepare the data
aqm.prepaffy = function(expressionset)
{
  pp1 = try(preprocess(expressionset))
  dataPLM = try(fitPLM(pp1, background = FALSE, normalize = FALSE))
  if(class(pp1)=='try-error' || class(dataPLM)=='try-error')
    warning("RLE and NUSE plots from the package 'affyPLM' cannot be produced for this data set.")

  Affyobject = list("dataPLM" = dataPLM)
  class(Affyobject) = "aqmobj.prepaffy"
  return(Affyobject)
}

## phrase book
phrase = list(
  fig = function(x) sprintf("The figure <!-- FIG --> shows the <i>%s</i> plot.", x),
  preproc = function(x) sprintf("%s values are computed from the preprocessed data (after background correction and quantile normalisation).", x),
  affyPLM = "This diagnostic plot is based on tools provided in the <i>affyPLM</i> package."
  )

##RNA digestion
aqm.rnadeg = function(expressionset)
  {
    sN = sampleNames(expressionset)
    numArrays = length(sN)
    maxcol = 8
    colors = brewer.pal(maxcol, "Dark2")
    acol = if(numArrays <= maxcol)
      colors[seq_len(numArrays)] else   
      c(colors, sample(colors, numArrays-length(colors), replace = TRUE))
    acol = addOpacity(acol, numArrays)
    
    rnaDeg = function() {
      plotAffyRNAdeg(AffyRNAdeg(expressionset, log.it = TRUE),
                     lwd = 1, col=acol)
      legend("topright", lty=1, lwd=2, col=acol, legend = paste(seq(along=acol)))
    }
    
    annotation = namedEmptyList(numArrays)
    for(i in seq(along=annotation))
      annotation[[i]] = list(title = sprintf("Array %d: %s", i, sN[i]),
                             linkedids = names(annotation)[i])

    legend = paste(phrase$fig("RNA digestion"),
                   phrase$preproc("Shown"),
       "In this plot each array is represented by a single line; move the mouse over the lines to see their corresponding sample names. The plot can be used to identify array(s) that have a slope very different from the others. This could indicate that the RNA used for that array has been handled differently from what was done for the other arrays.",
                   phrase$affyPLM)

    title = "RNA digestion plot"
    section = "Affymetrix specific plots"
    
    out = list("plot" = rnaDeg,
               "section" = section,
               "title" = title,
               "legend" = legend,
               "shape" = list("h" = 5.5, "w" =7),
               "svg" = list(annotation=annotation, getfun=SVGAnnotation::getMatplotSeries))
    
    class(out) = "aqmobj.rnadeg"
    return(out)
  }


##RLE
aqm.rle = function(expressionset, dataprep, affyproc, intgroup, subsample = 10000, ...)
  {   
    if (affyproc$dataPLM@model.description$R.model$which.parameter.types[3] == 1){
      medianchip <- apply(coefs(affyproc$dataPLM), 1, median)
      M <- sweep(coefs(affyproc$dataPLM),1,medianchip,FUN='-')
    } else {
      stop("It doesn't appear that a model with sample effects was used.")
    }

    prepaf = dataprep
    prepaf$dat = M
    rle = try(aqm.boxplot(expressionset, prepaf, intgroup, subsample = subsample, ...)) ## TO DO: Not sure if the score computation is correct though

    legend = paste(phrase$fig("Relative Log Expression (RLE)"),
      phrase$preproc("RLE"),
      "Arrays whose boxes have larger spread or are centred at M = 0 are potentially problematic.",
      phrase$affyPLM)
    
    title = "Relative Log Expression plot"
    
    section = "Affymetrix specific plots"

    rle2 = try(Mbox(affyproc$dataPLM, plot = FALSE))
    rlemed = rle2$stats[3,]
    rleout = which(abs(rlemed) > 0.1)

    out = list("plot" = rle$plot, "section" = section, "title" = title, "legend" = legend, "scores" = rlemed, "outliers" = rleout, shape = "square")
    class(out) = "aqmobj.rle"
    return(out)
  }

##NUSE
aqm.nuse = function(expressionset, dataprep, affyproc, intgroup, subsample = 10000, ...)
{ 
  ##bwplot for PLMset
  compute.nuse <- function(which){
    nuse <- apply(affyproc$dataPLM@weights[[1]][which,],2,sum)
    1/sqrt(nuse)
  }
  model <- affyproc$dataPLM@model.description$modelsettings$model
  if (affyproc$dataPLM@model.description$R.model$which.parameter.types[3] == 1 & affyproc$dataPLM@model.description$R.model$which.parameter.types[1] == 0 ){
    grp.rma.se1.median <- apply(se(affyproc$dataPLM), 1,median,na.rm=TRUE)
    grp.rma.rel.se1.mtx <- sweep(se(affyproc$dataPLM),1,grp.rma.se1.median,FUN='/')
  } else {
    which <-indexProbesProcessed(affyproc$dataPLM)
    ses <- matrix(0,length(which) ,4)
    if (affyproc$dataPLM@model.description$R.model$response.variable == 1){
      for (i in 1:length(which))
        ses[i,] <- compute.nuse(which[[i]])
    } else {
      stop("Sorry I can't currently impute NUSE values for this PLMset object")
    }
    grp.rma.se1.median <- apply(ses, 1,median)
    grp.rma.rel.se1.mtx <- sweep(ses,1,grp.rma.se1.median,FUN='/')
  }

  prepaf = dataprep
  prepaf$dat = grp.rma.rel.se1.mtx
  nuse = try(aqm.boxplot(expressionset, prepaf, intgroup, subsample = subsample, ...)) ## TO DO: Not sure if the score computation is correct though

  nuse$legend = paste(phrase$fig("Normalized Unscaled Standard Error (NUSE)"),
    phrase$preproc("NUSE"),
    "Potential low quality arrays are those whose box is substantially elevated or more spread out relative to the other arrays. NUSE values are comparable within one data set, but usually not across different data sets.",
    phrase$affyPLM)

  nuse$title = "Normalized Unscaled Standard Error (NUSE) plot"
  nuse$section = "Affymetrix specific plots"
  
  nuse$shape = "square"
  class(nuse) = "aqmobj.nuse"
  return(nuse) 
}


##QCStats
aqm.qcstats = function(expressionset, ...)
  {                
    qcStats = function() {
      plot.qc.stats(qc(expressionset, ...))
    }

    legend =  paste(phrase$fig("the Affymetrix recommended diagnostic"),
      "Please see the vignette of the package <i>simpleaffy</i> for a full explanation of the elements shown in this plot. Any metrics that is shown in red is outside the manufacturer's specified boundaries and suggests a potential problem; metrics shown in blue are considered acceptable.")

    title = "Diagnostic plot recommended by Affymetrix"
    section = "Affymetrix specific plots"

    out = list("plot" = qcStats, "section" = section, "title" = title, "legend" = legend, shape = "square")
    class(out) = "aqmobj.qcs"
    return(out)    
  }

##PM / MM
aqm.pmmm = function(expressionset, ...)
{  
  PM = density(as.matrix(log2(pm(expressionset))), ...)
  MM = density(as.matrix(log2(mm(expressionset))), ...)
  PMMM = list(PM = PM, MM = MM)

  legend = "Figure <!-- FIG --> shows the density distributions of the log<sub>2</sub> intensities grouped by the matching type of the probes. The blue line shows a density estimate (smoothed histogram) from intensities of perfect match probes (PM), the grey line, one from the mismatch probes (MM). We expect that MM probes have poorer hybridization than PM probes, and thus that the PM curve be to the right of the MM curve."

  title = "Perfect matches and mismatches"
  section = "Affymetrix specific plots"
  out = list("plot" = PMMM, "section" = section, "title" = title, "legend" = legend, shape = "square")
  class(out) = "aqmobj.pmmm"
  return(out)    
}



