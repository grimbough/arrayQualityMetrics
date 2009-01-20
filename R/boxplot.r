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
      box = bwplot(dataprep$dat ~ as.vector(col(dataprep$dat)), pch = "|", col = "black", do.out = FALSE, fill = "#1F78B4", horizontal = FALSE, box.ratio = 2, xlab = "", ylab = "Intensities", asp = "iso", ...)
      shape = "square"
    }
  legspe = if(dataprep$nchannels == 2) "The left panel corresponds to the red channel. The middle panel shows the green channel. The right panel shows the boxplots of log2(ratio)." else ""
  
  legend = sprintf("<b>Figure FIG</b> presents boxplots of the log<sub>2</sub>(Intensities). Each box corresponds to one array. %s Typically, one expects the boxes to have similar size (IQR) and y position (median).", legspe)

  title = "Boxplots"
  type = "Homogeneity between arrays"

  b = boxplot(dataprep$dat, plot=F, range=0)  
  bmeanstat = boxplot.stats(b$stat[3,])
  bmeanout = sapply(seq_len(length(bmeanstat$out)), function(x) which(b$stat[3,] == bmeanstat$out[x]))
  
  biqr = b$stat[4,] -  b$stat[2,] 
  biqrstat = boxplot.stats(biqr)
  biqrout = sapply(seq_len(length(biqrstat$out)), function(x) which(biqr == biqrstat$out[x]))
  
  out = list("plot" = box, "type" = type,  "title" = title, "legend" = legend, "scores" = list("mean" = b$stat[3,], "IQR" = biqr), "outliers" = list("mean" = bmeanout, "IQR" = biqrout), "shape" = shape)
  class(out) = "aqmobj.box"
  return(out)   
}
