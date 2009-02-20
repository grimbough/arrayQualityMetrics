aqm.boxplot = function(dataprep, ...)
{
  if(dataprep$nchannels == 2)
    {
      boxred = bwplot(dataprep$rc ~ as.vector(col(dataprep$rc)),pch = "|", col = "black", do.out = FALSE, fill = "#E31A1C", horizontal = FALSE, box.ratio = 2, xlab = "", ylab = "Red intensities", asp = "iso", ...)
      boxgreen = bwplot(dataprep$gc ~ as.vector(col(dataprep$gc)), pch = "|", col = "black", do.out = FALSE, fill = "#33A02C", horizontal = FALSE, box.ratio = 2,  xlab = "", ylab = "Green intensities", asp = "iso", ...)
      boxblue = bwplot(dataprep$dat ~ as.vector(col(dataprep$dat)), pch = "|", col = "black", do.out = FALSE, fill = "#1F78B4", horizontal = FALSE, box.ratio = 2,  xlab = "", ylab = "Log(Ratio)", asp = "iso", ...)
      box = list("boxred" = boxred, "boxgreen" = boxgreen, "boxblue" = boxblue)
      shape = "rect"
    } else  {
      box = bwplot(dataprep$dat ~ as.vector(col(dataprep$dat)), pch = "|", col = "black", do.out = FALSE, fill = "#1F78B4", horizontal = FALSE, box.ratio = 2, xlab = "", ylab = "Intensities", ...)
      shape = "square"
    }
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
