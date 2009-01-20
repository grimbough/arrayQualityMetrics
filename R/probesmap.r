aqm.probesmap = function(expressionset, dataprep, ...)
{  
  ht = featureData(expressionset)$hasTarget  

  df = data.frame("ProbesMap" = rep(ht,ncol(dataprep$dat)), "dat" = as.numeric(dataprep$dat))

  den = densityplot(~ dat, data = df, groups = ProbesMap, plot.points = FALSE, auto.key = list(lines=TRUE, points=FALSE, x=1, y=1, corner=c(1,1)), xlab = "Intensity")

  legend = "<b>Figure FIG</b> shows the density distributions of the log<sub>2</sub> ratios grouped by the mapping of the probes. Blue, density estimate of log<sub>2</sub> ratios of probes annotated \"TRUE\" in the <b>\"hasTarget\"</b> slot. Gray, probes annotated \"FALSE\" in the <b>\"hasTarget\"</b> slot."
  
  title = "Probes mapping"
  type = "Platform quality assessment"
  
  out = list("plot" = den, "type" = type, "title" = title, "legend" = legend, "shape" = "square")
  class(out) = "aqmobj.probesmap"
  return(out)   
}
