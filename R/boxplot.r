
within = function(x, range) {
  stopifnot(length(range)==2)
  return((x >= range[1]) & (x <= range[2]))
}

aqm.boxplot = function(expressionset, dataprep, intgroup = "Covariate", grouprep = FALSE, ...)
{

  b = boxplot(dataprep$dat, plot=FALSE, range=0)

  ## Outlier detection, applied to the boxplot statistics: to the medians (i.e. box positions),
  ## and differences between upper hinge and lower hinge (i.e. box sizes)
  bmedian = b$stat[3,]
  bmedianstat = boxplot.stats(bmedian)
  biqr = b$stat[4,] -  b$stat[2,] 
  biqrstat = boxplot.stats(biqr)

  bmedianout = !within(bmedian, bmedianstat$stats[c(1,5)])
  biqrout    = !within(biqr,    biqrstat$stats[c(1,5)])
  outliers   = bmedianout | biqrout

  
  bsdat   = boxplot(dataprep$dat, plot=FALSE)

  stopifnot( dataprep$nchannels %in% c(1,2) )
  if(dataprep$nchannels == 2)  {
    bsrc = boxplot(dataprep$rc, plot=FALSE)
    bsgc = boxplot(dataprep$gc, plot=FALSE)
    
    bigdat  = c(as.numeric(bsrc$stats), as.numeric(bsgc$stats), as.numeric(bsdat$stats))
    colbigd = c(as.vector(col(bsrc$stat)), as.vector(col(bsgc$stats)), as.vector(col(bsdat$stats)))
    fac = c(rep("a. Red Channel",   prod(dim(bsrc$stats ))),
            rep("b. Green Channel", prod(dim(bsgc$stats ))),
            rep("c. Log(Ratio)",    prod(dim(bsdat$stats))))
  } else {
    bigdat  = as.numeric(bsdat$stats)
    colbigd = as.vector(col(bsdat$stats))
  }
  
  if(all(intgroup %in% names(phenoData(expressionset)@data)) && grouprep) {
    gp = intgroup[1]
    gpcont = pData(expressionset)[colnames(pData(expressionset))==gp,]
    coloursb = brewer.pal(8, rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[6])
    intgroupcont = gpcont[,gp]
    colsb = coloursb[as.factor(intgroupcont)]
    colsb[outliers] = "grey"
     
    key = list(rect = list(col=unlist(coloursb)[as.factor(levels(as.factor(unlist(gpcont[,1]))))]),
               text = list(levels(as.factor(unlist(gpcont[,1])))))
    key$rep = FALSE
    foo = draw.key(key = key)
    
  } else {
    colsb = rep("#1F78B4", dataprep$numArrays) 
    colsb[outliers] = "grey"
    
    foo = NULL
  }
    
  if(dataprep$nchannels == 2) {
    formula = colbigd ~ bigdat | factor(fac)
    lay = c(3,1)
  } else {
    formula = colbigd ~ bigdat
    lay = c(1,1)
  }  

  box = bwplot(formula, layout=lay, as.table=TRUE,
        strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
        horizontal = TRUE, groups = colbigd, 
        pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
        xlab = "", fill = colsb, panel = panel.superpose, panel.groups = panel.bwplot, main=foo)


  shape = list("h" = if(dataprep$numArrays > 50) (dataprep$numArrays/8) else 10, "w" = 10) 
  
  legspe = if(dataprep$nchannels == 2) "Left panel: red channel, middle panel: green channel, right panel: log<sub>2</sub>(ratio)." else ""
  
  legend = sprintf("The figure <!-- FIG --> presents boxplots of the log<sub>2</sub>(Intensities). Each box corresponds to one array. %s It gives a simple summary of the distribution of probe intensities across all arrays. Typically, one expects the boxes to have similar size (IQR) and y position (median). If the distribution of an individual array is very different from the others, this may indicate an experimental problem. After normalisation, the distributions should be similar.", legspe)

  title = "Boxplots"
  section = "Array intensity distributions"

  out = list("plot" = box,
             "section" = section,
             "title" = title,
             "legend" = legend,
             "scores"   = list("median" = bmedian, "IQR" = biqr),
             "outliers" = list("median" = bmedianout, "IQR" = biqrout),
             "shape" = shape)
  class(out) = "aqmobj.box"
  
  return(out)   
}
