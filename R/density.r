dens = function(obj)
  {
    dlist = apply(obj, 2, function(x){
      rg = quantile(x, na.rm = TRUE, probs = c(0.01, 0.99))
      density(x, na.rm = TRUE, from = rg[1], to = rg[2])
    })
    names(dlist) = seq_along(dlist)
    ddf = do.call(make.groups, lapply(dlist, function(l) with(l, data.frame(x = x, y = y))))
    return(ddf)
  }

namedEmptyList = function(n) {
  x = vector(mode="list", length = n)
  names(x) = sprintf("aqm_%d", seq(along=x))
  return(x)
}
    
aqm.density = function(x)
{
  
  cl = intgroupColours(x)
  lwd = lty = 1
  name = "density"
  
  if(x$nchannels == 2)
    {  
      den1 = dens(x$R)
      den2 = dens(x$G)
      den3 = dens(x$M)
      ddf  = rbind(den1, den2, den3)
      panels = factor(rep(1:3, times = c(nrow(den1), nrow(den2), nrow(den3))),
        levels = 1:3,
        labels = c("a. Red Channel", "b. Green Channel", "c. Log2(Ratio)"))
      formula = y ~ x | panels
      lay = c(3,1)

      ## the mappings between report objects (arrays) and plot objects (lines)
      n = x$numArrays
      getPlotObjIdFromReportObjId = sprintf("
   // report object id -> plot object id
   function(r) {
     var a = parseInt(r.replace(/^r:/, ''));
     var n = %d;
     if(isNaN(a)) throw new Error('Invalid report object id ' + x);
     return(['p:' + a, 'p:' + (a+n), 'p:' + (a+2*n)]);
   }", n)

      getReportObjIdFromPlotObjId = function(x) {
        j = as.integer(sub("^p:", "", x))
        stopifnot(length(j)==1, !is.na(j), j>0)
        paste("r", (j-1) %% n + 1, sep = ":")
      }

      svgPar = new("svgParameters",
                   name           = name,
                   numPlotObjects = 3L*x$numArrays,
                   getPlotObjIdFromReportObjId = getPlotObjIdFromReportObjId,
                   getReportObjIdFromPlotObjId = getReportObjIdFromPlotObjId)

      
    } else {  ## nchannels==1
      ddf = dens(x$M)
      formula = y ~ x
      lay = c(1,1)
      svgPar = new("svgParameters",
                   name = name,
                   numPlotObjects = x$numArrays)
    }

  
  den = xyplot(formula, ddf, groups = which, layout = lay,
    type = "l", ylab = "Density", xlab="",
    main = if(!is.null(cl$key)) draw.key(key = cl$key),
    strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
    scales = list(relation="free"),
    col = cl$arrayColours, lwd = lwd, lty = lty)
  
  legend = "The figure <!-- FIG --> shows density estimates (smoothed histograms) of the data. Typically, the distributions of the arrays should have similar shapes and ranges. Arrays whose distributions are very different from the others should be considered for possible problems. Various features of the distributions can be indicative of quality related phenomena. For instance, high levels of background will shift an array's distribution to the right. Lack of signal diminishes its right right tail. A bulge at the upper end of the intensity range often indicates signal saturation."
  
  new("aqmReportModule",
      plot    = den,
      section = "Array intensity distributions",
      title   = "Density plots",
      legend  = legend,
      size    = c(w = 3+3*lay[1], h = 5*lay[2]),
      svg     = svgPar
      ) ## new
}
