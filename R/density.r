dens = function(obj)
  {
    dlist = apply(obj, 2, density, na.rm = TRUE)
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
  
  if(x$nchannels == 2)
    {  
      den1 = dens(x$rc)
      den2 = dens(x$gc)
      den3 = dens(x$dat)
      ddf  = rbind(den1, den2, den3)
      panels = factor(rep(1:3, each = c(nrow(den1), nrow(den2), nrow(den3))),
        levels = 1:3,
        labels = c("a. Red Channel", "b. Green Channel", "c. Log2(Ratio)"))
      formula = y ~ x | panels
      lay = c(3,1)

      ## the mappings between report objects (arrays) and plot objects (lines)
      idFun = sprintf(
"[ // report object id -> plot object id
   function(x) {
     var a = parseInt(x.replace('^ro:', ''));
     var n = %d;
     ['po:' + a, 'po:' + (a+n), 'po:' + (a+2*n)];
   },
   // plot object id -> report object id
   function(x) {
     var r = parseInt(x.replace('^po:', ''));
     'ro:' + (r-1) %% %d + 1;
   }
]", x$numArrays, x$numArrays)

      svgPar = if(x$usesvg)
        new("svgParameters",
              defined = TRUE,
              numPlotObjects = 3*x$numArrays,
              idFun = idFun) else new("svgParameters")
    } else {  ## nchannels==1
      ddf = dens(x$dat)
      formula = y ~ x
      lay = c(1,1)
      svgPar = if(x$usesvg)
        new("svgParameters",
              defined = TRUE,
              numPlotObjects = x$numArrays) else new("svgParameters")
    }

  
  den = xyplot(formula, ddf, groups = which, layout = lay,
    type = "l", ylab = "Density", xlab="",
    main = if(!is.null(cl$key)) draw.key(key = cl$key),
    strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
    scales = list(relation="free"),
    col = cl$arrayColours, lwd = lwd, lty = lty)
  
  legend = "The figure <!-- FIG --> shows density estimates (smoothed histograms) of the data. Typically, the distributions of the arrays should have similar shapes and ranges. Arrays whose distributions are very different from the others should be considered for possible problems. Move the mouse over the lines in the plot to see the  corresponding sample names.%s<BR>On raw data, a bimodal distribution can be indicative of an array containing a spatial artefact; a distribution shifted to the right of an array with abnormally high background intensities."
  
  new("aqmReportModule",
      "plot"    = den,
      "section" = "Array intensity distributions",
      "title"   = "Density plots",
      "legend"  = legend,
      "shape"   = list("h" = 5, "w" = 3+3*lay[1]),
      "svg"     = svgPar
      ) ## new
}
