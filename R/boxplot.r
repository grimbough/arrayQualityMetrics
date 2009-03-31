aqm.boxplot = function(expressionset, dataprep, intgroup = "Covariate", grouprep = FALSE, ...)
{

  if(dataprep$nchannels == 2)
    {
      bsrc = boxplot(dataprep$rc, plot=F)
      bsgc = boxplot(dataprep$gc, plot=F)
      bsdat = boxplot(dataprep$dat, plot=F)
          
      bigdat = c(as.numeric(bsrc$stats), as.numeric(bsgc$stats), as.numeric(bsdat$stats))
      colbigd = c(as.vector(col(bsrc$stat)), as.vector(col(bsgc$stats)), as.vector(col(bsdat$stats)))
      fac = c(rep("a. Red Channel",dim(bsrc$stats)[1]*dim(bsrc$stats)[2]), rep("b. Green Channel",dim(bsgc$stats)[1]*dim(bsgc$stats)[2]), rep("c. Log(Ratio)",dim(bsdat$stats)[1]*dim(bsdat$stats)[2]))

    }
  else {
    bsdat = boxplot(dataprep$dat, plot=F)
    bigdat = as.numeric(bsdat$stats)
    colbigd = as.vector(col(bsdat$stats))
  }

  
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
          igcbd = intgroupcont[colbigd]
          box = bwplot(bigdat ~ colbigd | factor(fac), horizontal=F, groups = igcbd, layout=c(3,1), as.table=T, strip = function(..., bg) strip.default(..., bg ="#cce6ff"), pch = "|",  col = "black", do.out = FALSE, box.ratio = 2, xlab = "", fill = colsb, panel = panel.superpose, panel.groups = panel.bwplot, main=foo)
        }
      else {
        igcbd = intgroupcont[colbigd]
        box = bwplot(bigdat ~ colbigd, horizontal=F, groups = igcbd, pch = "|",  col = "black", do.out = FALSE, box.ratio = 2, xlab = "", fill = colsb, panel = panel.superpose, panel.groups = panel.bwplot, main=foo, ylab = "")
      }
    } else {
      if(dataprep$nchannels == 2)
        {
         box = bwplot(bigdat ~ colbigd | factor(fac), groups = which, horizontal=F, layout=c(3,1), as.table=T, strip = function(..., bg) strip.default(..., bg ="#cce6ff"), pch = "|",  col = "black", do.out = FALSE, box.ratio = 2, xlab = "", fill = c("#E31A1C", "#33A02C", "#1F78B4"))
         ## Need to find a way to have the colours on the plot
         ##E31A1C red
         ##33A02C green
         ##1F78B4 blue
        } else box = bwplot(bigdat ~ colbigd, horizontal=F, pch = "|",  col = "black", do.out = FALSE, box.ratio = 2, xlab = "", ylab= "",  fill = "#1F78B4")
    }

  if(dataprep$nchannels == 2)
    shape = "rect" else shape = "square"
  
  legspe = if(dataprep$nchannels == 2) "The left panel corresponds to the red channel. The middle to the green channel. The right panel shows the boxplots of log<sub>2</sub>(ratio)." else ""
  
  legend = sprintf("The figure <!-- FIG --> presents boxplots of the log<sub>2</sub>(Intensities). Each box corresponds to one array. %s It gives a simple summary of the distribution of probe intensities across all arrays. Typically, one expects the boxes to have similar size (IQR) and y position (median). If the distribution of an individual array is very different from the others, this may indicate an experimental problem. After normalisation, the distributions should be similar.", legspe)

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
