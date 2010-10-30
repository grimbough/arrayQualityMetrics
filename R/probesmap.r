aqm.probesmap = function(x)
{

  if(!("hasTarget" %in% colnames(x$fData)))
    return(NULL)
  
  df = data.frame(
    "hasTarget" = rep(x$fData$hasTarget, ncol(x$dat)),
    "dat"       = as.vector(x$dat))

  den = densityplot( ~ dat,
    data = df,
    groups = hasTarget,
    plot.points = FALSE,
    auto.key = list(lines = TRUE, points = FALSE, x = 1, y = 0.99, corner=c(1,1)),
    xlab = "Intensity")
  
  legend = "The figure <!-- FIG --> shows the intensity distributions (density estimates, pooled across all arrays) grouped by the value of the feature-level property <tt>hasTarget</tt>. This plot allows you to see whether the intensities of different types of features are systematically different. For instance, if <tt>hasTarget</tt> indicates whether or not a feature maps to a known gene, then you would expect the distribution of those features for which this property is true to be shifted towards higher values."
    
  new("aqmReportModule",
      plot = den,
      section = "Feature stratification",
      title = "Feature stratification",
      legend = legend,
      shape = list("h"=6,"w"=6))
}
