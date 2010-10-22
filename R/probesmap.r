aqm.probesmap = function(x)
{

  if(!("hasTarget" %in% colnames(fData(x$expressionset))))
    return(NULL)
  
  if(inherits(x$expressionset, "BeadLevelList")) {
    warning("probe mapping plot is not implemented for BeadLevelList objects.")
    return(NULL)
  } 

  df = data.frame(
    "ProbesMap" = rep(featureData(x$expressionset)$hasTarget, ncol(x$dat)),
    "dat"       = x$dat)

  den = densityplot( ~ dat,
    data = df,
    groups = ProbesMap,
    plot.points = FALSE,
    auto.key = list(lines = TRUE, points = FALSE, x = 1, y = 0.99, corner=c(1,1)),
    xlab = "Intensity")
  
  legend = "The figure <!-- FIG --> shows the density distributions of the log<sub>2</sub> ratios grouped by the mapping of the probes, i.e. by whether they have <tt>TRUE</tt> or <tt>FALSE</tt> in the <tt>hasTarget</tt> slot. The goal of this plot is to check whether the intensities of the probes mapping to a target are shifted to higher values compared to the intensities of the probes without known target. Two curves superposed may mean a problem in the annotation of the features."
    
  new("aqmTrellis",
      plot = den,
      section = "Probe stratification",
      title = "Probes mapping",
      legend = legend,
      shape = list("h"=6,"w"=6))
}
