dens = function(x)
  {
    rg = quantile(x, na.rm = TRUE, probs = c(0.02, 0.98))
    rg = rg + diff(rg)*c(-1,1)*0.04
    x[ (x<rg[1]) | (x>rg[2]) ] = NA

    dlist = apply(x, 2, density, na.rm = TRUE)
    names(dlist) = seq_along(dlist)

    ddf = do.call(make.groups, lapply(dlist, function(l) with(l, data.frame(x = x, y = y))))
    return(ddf)
  }

namedEmptyList = function(n) {
  x = vector(mode="list", length = n)
  names(x) = sprintf("aqm_%d", seq(along=x))
  return(x)
}

aqm.density = function(x, ...)
{

  lwd = lty = 1

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

      getReportObjIdFromPlotObjId = function(j) {
        stopifnot(length(j)==1, !is.na(j), j>0L)
        (j-1L) %% n + 1L
      }

      svgPar = new("svgParameters",
                   numPlotObjects = 3L*x$numArrays,
                   getReportObjIdFromPlotObjId = getReportObjIdFromPlotObjId)


    } else {  ## nchannels==1
      ddf = dens(x$M)
      formula = y ~ x
      lay = c(1,1)
      svgPar = new("svgParameters",
                   numPlotObjects = x$numArrays)
    }


  den = xyplot(formula, ddf, groups = which, layout = lay,
    type = "l", ylab = "Density", xlab="",
    main = if(!is.null(x$key)) draw.key(key = x$key),
    strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
    scales = list(relation="free"),
    col = x$arrayColors, lwd = lwd, lty = lty)

  legend = "The figure <!-- FIG --> shows density estimates (smoothed histograms) of the data. Typically, the distributions of the arrays should have similar shapes and ranges. Arrays whose distributions are very different from the others should be considered for possible problems. Various features of the distributions can be indicative of quality related phenomena. For instance, high levels of background will shift an array's distribution to the right. Lack of signal diminishes its right right tail. A bulge at the upper end of the intensity range often indicates signal saturation."

  new("aqmReportModule",
      plot    = den,
      section = "Array intensity distributions",
      title   = "Density plots",
      id      = "dens",
      legend  = legend,
      size    = c(w = 2+2*lay[1], h = 3.5*lay[2] + length(x$key$rect$col) * 0.2),
      colors  = x$arrayColors,
      svg     = svgPar)
}
