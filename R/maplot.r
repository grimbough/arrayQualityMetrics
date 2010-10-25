##Function to perform the MAplot
aqm.maplot = function(x) {
  
  app = 4 + 2*(sum(x$numArrays>c(4,6)))
  nfig = ceiling(x$numArrays/8)
    
  xlimMA = quantile(x$A, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
  ylimMA = quantile(x$M, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)

  i = seq_len(x$numArrays)
  dummy.df = data.frame(
    i = factor(i, levels = i),
    x = i,
    y = i)
  
  trobj = xyplot(y ~ x | i, dummy.df,
    xlim = xlimMA,
    ylim = ylimMA,
    xlab = "A",
    ylab = "M",
    panel = function(x, y, ...) panel.smoothScatter(x=x$A[, x], y=x$M[, y], nbin = 250, raster=TRUE, ...),
    as.table = TRUE,      
    layout = c(app/2, 2, 1),
    asp = "iso",
    strip = function(..., bg) strip.default(..., bg ="#cce6ff")
    )
  
  id.firstpage = seq_len(app)
  
  ma = lapply(seq_len(nfig), function(i)
    {
      id.thispage = (i-1) * app + id.firstpage
      id.thispage = id.thispage[id.thispage <= x$numArrays]
      update(trobj, index.cond = list(id.thispage))
    })
  
  legspe1 = if(x$nchannels == 1) "where I<sub>1</sub> is the intensity of the array studied, and I<sub>2</sub> is the intensity of a \"pseudo\"-array that consists of the median across arrays." else "where I<sub>1</sub> and I<sub>2</sub> are the intensities of the two channels."
    
  legspe2 = if(x$classori == "BeadLevelList") "The calculations are done on the summarized data obtained by using the function createBeadSummaryData from the package beadarray." else ""
    
  legend = sprintf("The figure <!-- FIG --> shows the MA plot for each array. M and A are defined as :<br>M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>)),<br>%s %s Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, and there should be no trend in M as a function of A. If there is a trend in the lower range of A, this often indicates that the arrays have different background intensities; this may be addressed by background correction. A trend in the upper range of A can indicate saturation of the measurements; in mild cases, this may be addressed by non-linear normalisation (e.g. quantile normalisation).", legspe1, legspe2)

  mamean = colMeans(abs(x$M), na.rm=TRUE)
  mastat = boxplot.stats(mamean)
  maout  = sapply(seq_len(length(mastat$out)), function(x) which(mamean == mastat$out[x])) ## TODO this can be done more elegantly
  
  new("aqmReportModule",
      plot = ma,
      section = "Individual array quality",
      title = "MA plots",
      legend = legend,
      outliers = maout,
      shape = list(h=6, w=10))
}

