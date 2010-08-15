dens = function(obj, ...)
  {
    dlist <- apply(obj, 2, density, na.rm = TRUE)
    names(dlist) <- seq_along(dlist)
    ddf <- do.call(make.groups, lapply(dlist, function(l) with(l, data.frame(x = x, y = y))))
    return(ddf)
  }

aqm.density = function(expressionset, dataprep, intgroup, outliers = c(), ...)
{
  sN = sampleNames(expressionset)

  if(dataprep$nchannels == 2) {  
    den1 = dens(dataprep$rc)
    den2 = dens(dataprep$gc)
    den3 = dens(dataprep$dat)
    ddf  = rbind(den1, den2, den3)
    panels = factor(rep(1:3, each = c(nrow(den1), nrow(den2), nrow(den3))),
      levels = 1:3,
      labels = c("a. Red Channel", "b. Green Channel", "c. Log2(Ratio)"))
    formula = y ~ x | panels
    lay = c(3,1)

    annotation = vector(mode="list", length = 3*length(sN))
    names(annotation) = sprintf("line%d", seq(along=annotation))
    for(i in seq(along=annotation)) {
      s = ((i-1) %% length(sN))+1
      annotation[[i]] = list(title=sN[s], linkedids=names(annotation)[s+(0:2)*length(sN)])
    }

  } else {
    ddf = dens(dataprep$dat)
    formula = y ~ x
    lay = c(1,1)

    annotation = vector(mode="list", length = length(sN))
    names(annotation) = sprintf("line%d", seq(along=annotation))
    for(i in seq(along=annotation))
      annotation[[i]] = list(title=sN[i], linkedids=names(annotation)[i])
  }
  
  cl = intgroupColours(intgroup, expressionset)

  transparencysuffix = "80"  ## FIXME
  if(length(outliers)>0)
    transparencysuffix[outliers] = "00"

  lwd = rep(1,dataprep$numArrays)
  lty = rep(1,dataprep$numArrays)
  
  den = xyplot(formula, ddf, groups = which, layout = lay,
    type = "l", ylab = "Density", xlab="",
    strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
    scales = list(relation="free"),
    col = cl$arrayColours, lwd = lwd, lty = lty, ...)
  
  shape = list("h" = 5, "w" = 2+lay[1]*4)
      
  legend = "The figure <!-- FIG --> shows density estimates (smoothed histograms) of the data. Typically, the distributions of the arrays should have similar shapes and ranges. Arrays whose distributions are very different from the others should be considered for possible problems. On raw data, a bimodal distribution can be indicative of an array containing a spatial artefact; a distribution shifted to the right of an array with abnormal higher background intensities."
  
  title = "Density plots"
  section = "Array intensity distributions"
  
  out = list("plot" = den,
             "section" = section,
             "title" = title,
             "legend" = legend,
             "shape" = shape,
             "svg" = list(annotation=annotation))
    
  class(out) = "aqmobj.dens"
  return(out)   
}
