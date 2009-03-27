mybwplot = function(a, intgroupcont = NULL, ylab, colsb, ...)
  {
    if(!is.null(intgroupcont)) bwplot(a ~ as.vector(col(a)), groups = intgroupcont[as.vector(col(a))], pch = "|",  col = "black", do.out = FALSE,  horizontal = FALSE, box.ratio = 2, xlab = "", ylab = ylab, asp = "iso", fill = colsb, panel = panel.superpose, panel.groups = panel.bwplot, ...) else bwplot(a ~ as.vector(col(a)), pch = "|",  col = "black", do.out = FALSE,  horizontal = FALSE, box.ratio = 2, xlab = "", ylab = ylab, asp = "iso", fill = colsb, ...)
  }

aqm.boxplot = function(expressionset, dataprep, intgroup = "Covariate", grouprep = FALSE, ...)
{
  if(all(intgroup %in% names(phenoData(expressionset)@data)) && grouprep == TRUE)
      {

      gp = intgroup[1]
      gpcont = pData(expressionset)[colnames(pData(expressionset))==gp]
      coloursb = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[6])
      intgroupcont = gpcont[,gp]
      colsb = coloursb[as.factor(intgroupcont)]
     
      key = list(rect = list(col=colsb[as.factor(levels(as.factor(unlist(intgroupcont))))]), text = list(levels(as.factor(unlist(intgroupcont)))))
      key$rep = FALSE
       
      foo = draw.key(key = key)

      if(dataprep$nchannels == 2)
        {
          boxred = mybwplot(dataprep$rc, intgroupcont = intgroupcont, ylab = "Red Intensities", colsb = colsb, main = foo, ...)
          boxgreen = mybwplot(dataprep$gc, intgroupcont = intgroupcont, ylab = "Green intensities", colsb = colsb, main = foo, ...)
          boxblue = mybwplot(dataprep$dat, intgroupcont = intgroupcont, ylab = "Log(Ratio)", colsb = colsb, main = foo, ...) 
        }
      else box = mybwplot(dataprep$dat, intgroupcont = intgroupcont, ylab = "Intensities", colsb = colsb, main = foo, ...)
    } else {
      if(dataprep$nchannels == 2)
        {
          boxred = mybwplot(dataprep$rc, ylab = "Red intensities", colsb = "#E31A1C", ...) 
          boxgreen = mybwplot(dataprep$gc, ylab = "Green intensities", colsb = "#33A02C", ...) 
          boxblue = mybwplot(dataprep$dat, ylab = "Log(Ratio)",  colsb = "#1F78B4", ...) 
        } else box = mybwplot(dataprep$dat, ylab = "Intensities",  colsb = "#1F78B4", ...) 
    }

  if(dataprep$nchannels == 2)
    {
      box = list("boxred" = boxred, "boxgreen" = boxgreen, "boxblue" = boxblue)
      shape = "rect"} else shape = "square"
  
  legspe = if(dataprep$nchannels == 2) "The left panel corresponds to the red channel. The middle to the green channel. The right panel shows the boxplots of log<sub>2</sub>(ratio)." else ""
  
  legend = sprintf("The figure <!-- FIG --> presents boxplots of the log<sub>2</sub>(Intensities). Each box corresponds to one array. %s It gives a simple summary of the distribution of probe intensities accross all arrays. Typically, one expects the boxes to have similar size (IQR) and y position (median). If the distribution of an individual array is very different from the others, this may indicate an experimental problem. After normalisation, the distributions should be similar.", legspe)

  title = "Boxplots"
  section = "Array intensity distributions"

  b = boxplot(dataprep$dat, plot=F, range=0)  
  bmeanstat = boxplot.stats(b$stat[3,])
  bmeanout = sapply(seq_len(length(bmeanstat$out)), function(x) which(b$stat[3,] == bmeanstat$out[x]))
  
  biqr = b$stat[4,] -  b$stat[2,] 
  biqrstat = boxplot.stats(biqr)
  biqrout = sapply(seq_len(length(biqrstat$out)), function(x) which(biqr == biqrstat$out[x]))
  
  out = list("plot" = box, "section" = section, "title" = title, "legend" = legend, "scores" = list("mean" = b$stat[3,], "IQR" = biqr), "outliers" = list("mean" = bmeanout, "IQR" = biqrout), "shape" = shape)
  class(out) = "aqmobj.box"
  return(out)   
}
