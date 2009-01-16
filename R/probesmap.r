aqm.probesmap = function(expressionset, dataprep, ...)
{  
  key = list(lines = list(col=c("blue","grey")), text = list(c("Mapped","Unmapped")), space="right")
  
  ht = featureData(expressionset)$hasTarget  

  df = data.frame("ProbesMap" = rep(ht,ncol(dataprep$dat)), "dat" = as.numeric(dataprep$dat))

  den = densityplot(~ dat, data = df, groups = ProbesMap, plot.points = FALSE, key = key, col = c("grey","blue"), xlab = "Intensity", ...)

  legend = "<b>Figure FIG</b> shows the density distributions of the log<sub>2</sub> ratios grouped by the mapping of the probes. Blue, density estimate of log<sub>2</sub> ratios of probes annotated \"TRUE\" in the <b>\"hasTarget\"</b> slot. Gray, probes annotated \"FALSE\" in the <b>\"hasTarget\"</b> slot."
  
  title = "Probes mapping"
  type = "Platform quality assessment"
  
  out = list("plot" = den, "type" = type, "title" = title, "legend" = legend)
  class(out) = "aqmobj.probesmap"
  return(out)   
}
