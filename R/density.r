dens = function(obj, ...)
  {
    dlist <- apply(obj, 2, density)
    names(dlist) <- seq_along(dlist)
    ddf <- do.call(make.groups, lapply(dlist, function(l) with(l, data.frame(x = x, y = y))))
    densout = xyplot(y ~ x, ddf, groups = which, type = "l", ylab = "Density", ...)
    return(densout)
  }

aqm.density = function(dataprep, ...)
{  
  if(dataprep$nchannels == 2)
    {  
      den1 = dens(dataprep$rc, xlab = "Red Intensity")
      den2 = dens(dataprep$gc, xlab = "Green Intensity")
      den3 = dens(dataprep$dat, xlab = "Log(Ratio)", auto.key = list(lines=TRUE, points=FALSE, x=1, y=0.99, corner=c(1,1)))
      den = list("den1" = den1, "den2" = den2, "den3" = den3)
      shape = "rect"
    } else {
      den = dens(dataprep$dat, xlab = "Intensity", auto.key = list(lines=TRUE, points=FALSE, x=1, y=0.99, corner=c(1,1)))      
      shape = "square"
    }
      
  legend = "The figure <!-- FIG --> shows density estimates (smoothed histograms) of the data. Typically, the distributions of the arrays should have similar shapes and ranges. Arrays whose distributions are very different from the others should be considered for possible problems. On raw data, a bimodal distribution can be indicative of an array containing a spatial artifact and an array shifted to the right of an array with abnormal higher background intensities."
  
  title = "Density plots"
  section = "Array intensities distributions"
  
  out = list("plot" = den, "section" = section, "title" = title, "legend" = legend, "shape" = shape)
  class(out) = "aqmobj.dens"
  return(out)   
}
