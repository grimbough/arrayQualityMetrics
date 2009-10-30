dens = function(obj, ...)
  {
    dlist <- apply(obj, 2, density, na.rm = TRUE)
    names(dlist) <- seq_along(dlist)
    ddf <- do.call(make.groups, lapply(dlist, function(l) with(l, data.frame(x = x, y = y))))
    return(ddf)
  }

aqm.density = function(expressionset, dataprep, intgroup = "Covariate", grouprep = FALSE, ...)
{ 

  if(dataprep$nchannels == 2)
    {  

      den1 = dens(dataprep$rc)
      den2 = dens(dataprep$gc)
      den3 = dens(dataprep$dat)
      ddf = rbind(den1,den2,den3)
      fac = rep(c("a. Red Channel","b. Green Channel","c. Log(Ratio)"), each=nrow(den3))
    } else{
      ddf = dens(dataprep$dat)
    }

  if(all(intgroup %in% names(phenoData(expressionset)@data)) && grouprep == TRUE)
    {
      gp = intgroup[1]
      gpcont = pData(expressionset)[colnames(pData(expressionset))==gp]
      coloursb = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[6])
      intgroupcont = gpcont[,gp]
      colsb = coloursb[as.factor(intgroupcont)]
     
      key = list(rect = list(col=unlist(coloursb)[as.factor(levels(as.factor(unlist(gpcont[,1]))))]), text = list(levels(as.factor(unlist(gpcont[,1])))))
      key$rep = FALSE
  
      foo = draw.key(key = key)

      if(dataprep$nchannels == 2)
        {  
          den = xyplot(y ~ x | factor(fac), ddf, groups = rep(which,3), type = "l", ylab = "Density", xlab="", layout=c(3,1), strip = function(..., bg) strip.default(..., bg ="#cce6ff"), col = colsb, main=foo) #
          shape = "rect"
        } else {
          den = xyplot(y ~ x, ddf, groups = which, type = "l", ylab = "Density", xlab="", col = colsb, main=foo)
          shape = "square"
        }
    } else {  
      if(dataprep$nchannels == 2)
        {  
          den = xyplot(y ~ x | factor(fac), ddf, groups = rep(which,3), type = "l", ylab = "Density", xlab="", layout=c(3,1), strip = function(..., bg) strip.default(..., bg ="#cce6ff"), auto.key = list(lines=TRUE, points=FALSE, x=1, y=0.96, corner=c(1,1))) 
          shape = "rect"
        } else {
          den = xyplot(y ~ x, ddf, groups = which, type = "l", ylab = "Density", xlab="", auto.key = list(lines=TRUE, points=FALSE, x=1, y=0.99, corner=c(1,1)))      
          shape = "square"
        }
    }
      
  legend = "The figure <!-- FIG --> shows density estimates (smoothed histograms) of the data. Typically, the distributions of the arrays should have similar shapes and ranges. Arrays whose distributions are very different from the others should be considered for possible problems. On raw data, a bimodal distribution can be indicative of an array containing a spatial artefact and an array shifted to the right of an array with abnormal higher background intensities."
  
  title = "Density plots"
  section = "Array intensity distributions"
  
  out = list("plot" = den, "section" = section, "title" = title, "legend" = legend, "shape" = shape)
  class(out) = "aqmobj.dens"
  return(out)   
}
