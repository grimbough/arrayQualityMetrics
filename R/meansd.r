aqm.meansd = function(x, ...)
{
  msdplo = function() {
    meanSdPlot(x$M,
               cex.axis = 0.9,
               ylab = "Standard deviation of the intensities",
               xlab="Rank(mean of intensities)")
  }

  legend = "The figure <!-- FIG --> shows a density plot of the standard deviation of the intensities across arrays on the <i>y</i>-axis versus the rank of their mean on the <i>x</i>-axis. The red dots, connected by lines, show the running median of the standard deviation. After normalisation and transformation to a logarithm(-like) scale, one typically expects the red line to be approximately horizontal, that is, show no substantial trend. In some cases, a hump on the right hand of the x-axis can be observed and is symptomatic of a saturation of the intensities."

  new("aqmReportModule",
      plot    = msdplo,
      section = "Variance mean dependence",
      title   = "Standard deviation versus rank of the mean",
      id      = "msd",
      legend  = legend,
      colors  = x$arrayColors,
      size    = c(w=6, h=6))
}
