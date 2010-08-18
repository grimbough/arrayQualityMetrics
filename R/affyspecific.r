##Methods to prepare the data
aqm.prepaffy = function(expressionset, sN)
{
  pp1 = try(preprocess(expressionset))
  dataPLM = try(fitPLM(pp1, background = FALSE, normalize = FALSE))
  if(class(pp1)=='try-error' || class(dataPLM)=='try-error')
    warning("RLE and NUSE plots from the package 'affyPLM' cannot be produced for this data set.")

  Affyobject = list("dataPLM" = dataPLM, "sN" = sN)
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
    
    rnaDeg = function() {
      plotAffyRNAdeg(AffyRNAdeg(expressionset, log.it = TRUE),
                     lwd = 1, col=acol)
      legend("topright", lty=1, lwd=2, col=acol, legend = paste(seq(along=acol)))
    }
    
    annotation = vector(mode="list", length = numArrays)
    names(annotation) = sprintf("line%d", seq(along=annotation))
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
aqm.rle = function(affyproc, ...)
  {   
    rle = try(Mbox(affyproc$dataPLM, ylim = c(-1, 1),
      names = affyproc$sN, col = brewer.pal(9, "Set1")[2],
      whisklty = 0, staplelty = 0,
      main = "RLE", las = 3, cex.axis = 0.8, plot = FALSE, ...))

      legend = paste(phrase$fig("Relative Log Expression (RLE)"),
                   phrase$preproc("RLE"),
        "Arrays whose boxes have larger spread or are centred at M = 0 are potentially problematic.",
                   phrase$affyPLM)

    title = "Relative Log Expression plot"

    section = "Affymetrix specific plots"
    
    rlemed = rle$stats[3,]
    rleout = which(abs(rlemed) > 0.1)

    out = list("plot" = rle, "section" = section, "title" = title, "legend" = legend, "scores" = rlemed, "outliers" = rleout, shape = "square")
    class(out) = "aqmobj.rle"
    return(out)
  }

##NUSE
aqm.nuse = function(affyproc, ...)
{ 
  nuse = try(boxplot(affyproc$dataPLM, ylim = c(0.95, 1.5), names = affyproc$sN,
    outline = FALSE, col =  brewer.pal(9, "Set1")[2], main = "NUSE",
    las = 2, cex.axis = 0.8, plot=FALSE,...))

  legend = paste(phrase$fig("Normalized Unscaled Standard Error (NUSE)"),
                 phrase$preproc("NUSE"),
    "Potential low quality arrays are those whose box is substantially elevated or more spread out relative to the other arrays. NUSE values are comparable within one data set, but usually not across different data sets.",
                 phrase$affyPLM)

  title = "Normalized Unscaled Standard Error (NUSE) plot"
  section = "Affymetrix specific plots"

  nusemean = nuse$stat[3,]
  nusemeanstat = boxplot.stats(nusemean)
  nusemeanout = sapply(seq_len(length(nusemeanstat$out)), function(x) which(nuse$stat[3,] == nusemeanstat$out[x]))
  
  nuseiqr = nuse$stat[4,] -  nuse$stat[2,] 
  nuseiqrstat = boxplot.stats(nuseiqr)
  nuseiqrout = sapply(seq_len(length(nuseiqrstat$out)), function(x) which(nuseiqr == nuseiqrstat$out[x]))
  
  out = list("plot" = nuse, "section" = section, "title" = title, "legend" = legend, "scores" = list(mean = nusemean, iqr = nuseiqr), "outliers" = list(mean = nusemeanout, iqr = nuseiqrout), shape = "square")
  class(out) = "aqmobj.nuse"
  return(out) 
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


