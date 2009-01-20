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
      den3 = dens(dataprep$dat, xlab = "Log(Ratio)", auto.key = list(lines=TRUE, points=FALSE, x=1, y=1, corner=c(1,1)))
      den = list("den1" = den1, "den2" = den2, "den3" = den3)
      shape = "rect"
    } else {
      den = dens(dataprep$dat, xlab = "Intensity", auto.key = list(lines=TRUE, points=FALSE, x=1, y=1, corner=c(1,1)))      
      shape = "square"
    }
      
  legend = "<b>Figure FIG</b> shows density estimates (histograms) of the data. Arrays whose distributions are very different from the others should be considered for possible problems."
  
  title = "Density plots"
  type = "Homogeneity between arrays"
  
  out = list("plot" = den, "type" = type, "title" = title, "legend" = legend, "shape" = shape)
  class(out) = "aqmobj.dens"
  return(out)   
}
