aqm.probesmap = function(expressionset, dataprep, ...)
{  
  ht = featureData(expressionset)$hasTarget  

  df = data.frame("ProbesMap" = rep(ht,ncol(dataprep$dat)), "dat" = as.numeric(dataprep$dat))

  den = densityplot(~ dat, data = df, groups = ProbesMap, plot.points = FALSE, auto.key = list(lines=TRUE, points=FALSE, x=1, y=0.99, corner=c(1,1)), xlab = "Intensity")

  legend = "The <b>figure <!-- FIG --></b> shows the density distributions of the log<sub>2</sub> ratios grouped by the mapping of the probes. Blue, density estimate of log<sub>2</sub> ratios of probes annotated \"TRUE\" in the <b>\"hasTarget\"</b> slot. Gray, probes annotated \"FALSE\" in the <b>\"hasTarget\"</b> slot. The goal of this metric is to check whether the density of the intensities of the probes mapping for a target is shifted to the right compare to the density of intensities of the probes without known target. We could indeed expect that probes with a target will have, in overall, higher signal intensities. Two curves superposed may mean a problem in the annotation of the features. This metric can be used to compare two groups of probes in general, not necessarily the fact that they map for a target or not."
  
  title = "Probes mapping"
  type = "Probes stratification"
  
  out = list("plot" = den, "type" = type, "title" = title, "legend" = legend, "shape" = "square")
  class(out) = "aqmobj.probesmap"
  return(out)   
}
