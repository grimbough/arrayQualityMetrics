aqm.probesmap = function(x)
{

  if(!("hasTarget" %in% colnames(x$fData)))
    return(NULL)
  
  df = data.frame(
    "hasTarget" = rep(x$fData$hasTarget, ncol(x$A)),
    "A"         = as.vector(x$A))

  den = densityplot( ~ A,
    data = df,
    groups = hasTarget,
    plot.points = FALSE,
    auto.key = list(lines = TRUE, points = FALSE, x = 1, y = 0.99, corner=c(1,1)),
    xlab = "Intensity",
    xlim = quantile(df$A, probs=c(0.01, 0.99)), adjust=0.5)
  
  legend = "The figure <!-- FIG --> shows the intensity distributions (density estimates, pooled across all arrays) grouped by the value of the feature-level property <tt>hasTarget</tt>. This plot allows you to see whether the intensities of different types of features are systematically different. For instance, if <tt>hasTarget</tt> indicates whether or not a feature maps to a known gene, then you would expect the distribution of those features for which this property is <tt>TRUE</tt> to be shifted towards higher values."
    
  new("aqmReportModule",
      plot    = den,
      section = "Feature stratification",
      title   = "Feature stratification",
      id      = "feats",
      legend  = legend,
      colors  = x$arrayColors,
      size    = c(w=5, h=5))
}
