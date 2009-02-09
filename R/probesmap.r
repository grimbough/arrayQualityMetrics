aqm.probesmap = function(expressionset, dataprep, ...)
{  
  ht = featureData(expressionset)$hasTarget  

  df = data.frame("ProbesMap" = rep(ht,ncol(dataprep$dat)), "dat" = as.numeric(dataprep$dat))

  den = densityplot(~ dat, data = df, groups = ProbesMap, plot.points = FALSE, auto.key = list(lines=TRUE, points=FALSE, x=1, y=0.99, corner=c(1,1)), xlab = "Intensity")

  legend = "The figure <!-- FIG --> shows the density distributions of the log<sub>2</sub> ratios grouped by the mapping of the probes. Pink, density estimate of log<sub>2</sub> ratios of probes annotated \"TRUE\" in the <b>\"hasTarget\"</b> slot. Blue, probes annotated \"FALSE\" in the <b>\"hasTarget\"</b> slot. The goal of this metric is to check whether the density of the intensities of the probes mapping to a target is shifted to the right (meaning overall higher signal intensities) compare to the density of intensities of the probes without known target. Two curves superposed may mean a problem in the annotation of the features. This metric can be used to compare two groups of probes in general, not necessarily the fact that they map to a target or not."
  
  title = "Probes mapping"
  section = "Probes stratification"
  
  out = list("plot" = den, "section" = section, "title" = title, "legend" = legend, "shape" = "square")
  class(out) = "aqmobj.probesmap"
  return(out)   
}
