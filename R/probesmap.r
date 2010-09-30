aqm.probesmap = function(expressionset, dataprep, ...)
{

  if(!("hasTarget" %in% colnames(fData(expressionset))))
    return(NULL)
  
  if(inherits(expressionset, "BeadLevelList")) {
    warning("probe mapping plot is not implemented for BeadLevelList objects.")
    return(NULL)
  } 

  ht = featureData(expressionset)$hasTarget  

  df = data.frame(
    "ProbesMap" = rep(ht, ncol(dataprep$dat)),
    "dat" = as.numeric(dataprep$dat))

  den = densityplot( ~ dat,
    data = df,
    groups =  ProbesMap,
    plot.points = FALSE,
    auto.key = list(lines = TRUE, points = FALSE, x = 1, y = 0.99, corner=c(1,1)),
    xlab = "Intensity")
  
  legend = "The figure <!-- FIG --> shows the density distributions of the log<sub>2</sub> ratios grouped by the mapping of the probes, i.e. by whether they have <tt>TRUE</tt> or <tt>FALSE</tt> in the <tt>hasTarget</tt> slot. The goal of this plot is to check whether the intensities of the probes mapping to a target are shifted to higher values compared to the intensities of the probes without known target. Two curves superposed may mean a problem in the annotation of the features."
  
  title = "Probes mapping"
  section = "Probe stratification"
  
  out = list("plot" = den, "section" = section, "title" = title, "legend" = legend, "shape" = list("h"=6,"w"=6))
  class(out) = "aqmobj.probesmap"
  return(out)   
}
