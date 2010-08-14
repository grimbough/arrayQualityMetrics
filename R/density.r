dens = function(obj, ...)
  {
    dlist <- apply(obj, 2, density, na.rm = TRUE)
    names(dlist) <- seq_along(dlist)
    ddf <- do.call(make.groups, lapply(dlist, function(l) with(l, data.frame(x = x, y = y))))
    return(ddf)
  }

aqm.density = function(expressionset, dataprep, intgroup, outliers = c(), ...)
{ 

  if(dataprep$nchannels == 2) {  
    den1 = dens(dataprep$rc)
    den2 = dens(dataprep$gc)
    den3 = dens(dataprep$dat)
    ddf  = rbind(den1, den2, den3)
    panels = factor(rep(1:3, each = c(nrow(den1), nrow(den2), nrow(den3))),
      levels = 1:3,
      labels = c("a. Red Channel", "b. Green Channel", "c. Log2(Ratio)"))
    formula = sample_id ~ values | panels
    lay = c(3,1)
  } else {
    ddf = dens(dataprep$dat)
    formula = sample_id ~ values
    lay = c(1,1)
  }
  
  if (!missing(intgroup)) {
      gpcont = pData(expressionset)[colnames(pData(expressionset))==intgroup]
      coloursb = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[6])
      intgroupcont = gpcont[,intgroup]
      colsb = coloursb[as.factor(intgroupcont)]
      lwd = rep(1,dataprep$numArrays)
      lty = rep(1,dataprep$numArrays)
      if(length(outliers) != 0)
      {
         lwd[outliers] = 3
         colsb[-outliers]="grey"
         key2 = list(lines=list(col=colsb[outliers]), text = list(as.character(dataprep$sN)[outliers]), points=FALSE, corner=c(1,1), x=1, y =0.94)
         key2$rep = FALSE
      } else key2 = NULL
     
      key = list(rect = list(col=unlist(coloursb)[as.factor(levels(as.factor(unlist(gpcont[,1]))))]), text = list(levels(as.factor(unlist(gpcont[,1])))))
      key$rep = FALSE
  
      foo = draw.key(key = key)

      if(dataprep$nchannels == 2)
        {  
          den = xyplot(y ~ x | factor(fac), ddf, groups = rep(which,3), type = "l", ylab = "Density", xlab="", layout=c(3,1), strip = function(..., bg) strip.default(..., bg ="#cce6ff"), col = colsb, main=foo, lwd = lwd, lty = lty, key=key2, ...) 
          shape = "rect"
        } else {
          den = xyplot(y ~ x, ddf, groups = which, type = "l", ylab = "Density", xlab="", col = colsb, main=foo, lwd = lwd, lty = lty, ...)
          shape = "square"
        }
    } else {  
          lwd = rep(1,dataprep$numArrays)
	      lty = rep(1,dataprep$numArrays)
      if(length(outliers) != 0)
	  {
	  lwd[outliers] = 3
	  colo = rep("grey", dataprep$numArrays)
	  colo[outliers] = brewer.pal(8,"Set1")[seq_len(length(outliers))]	  
	  key = list(lines=list(col=colo[outliers]), text = list(as.character(dataprep$sN)[outliers]), points=FALSE, corner=c(1,1), x=1, y =0.94)
      key$rep = FALSE
      } else {
      colo = rep(brewer.pal(8,"Set1"), signif(dataprep$numArrays/8,0))[seq_len(dataprep$numArrays)]
	  key = list(lines=list(col=colo), text = list(as.character(dataprep$sN)), points=FALSE, corner=c(1,1), x=1, y =0.94)
      key$rep = FALSE
      }
      
      if(dataprep$nchannels == 2)
        {  
          den = xyplot(y ~ x | factor(fac), ddf, groups = rep(which,3), type = "l", ylab = "Density", xlab="", layout=c(3,1), strip = function(..., bg) strip.default(..., bg ="#cce6ff"),  lwd=lwd, lty=lty, col=colo, key=key, ...) 
          shape = "rect"
        } else {
          den = xyplot(y ~ x, ddf, groups = which, type = "l", ylab = "Density", xlab="", lwd=lwd, lty=lty, col=colo, key=key, ...) 
          shape = "square"
        }
    }
      
  legend = "The figure <!-- FIG --> shows density estimates (smoothed histograms) of the data. Typically, the distributions of the arrays should have similar shapes and ranges. Arrays whose distributions are very different from the others should be considered for possible problems. On raw data, a bimodal distribution can be indicative of an array containing a spatial artefact; a distribution shifted to the right of an array with abnormal higher background intensities."
  
  title = "Density plots"
  section = "Array intensity distributions"
  
  out = list("plot" = den, "section" = section, "title" = title, "legend" = legend, "shape" = shape)
  class(out) = "aqmobj.dens"
  return(out)   
}
