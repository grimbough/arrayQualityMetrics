##Function to perform the MAplot
maplotdraw = function(M, A, sN, numArrays, nchannels, class, ...)
  {
    app = 4 + 2*(sum(numArrays>c(4,6)))
    nfig = ceiling(numArrays/8)
    
    xlimMA = quantile(A, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
    ylimMA = quantile(M, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
    
    dummy.df = data.frame(sN = factor(sN, levels = sN),
      x = seq(along = sN),
      y = seq(along = sN))
    
    trobj = xyplot(y ~ x | sN, dummy.df,
      xlim = xlimMA,
      ylim = ylimMA,
      xlab = "A",
      ylab = "M",
      
      panel = function(x, y, ...) {
        x <- A[, x]
        y <- M[, y]
        ## FIXME: raster=TRUE would be nice to reduce the size of the currently enormous PDF files
        ##   but when tried last time (16.8.2010,  r52737), this made R crash  with
        ## *** caught segfault *** address 0x28, cause 'memory not mapped'...
        panel.smoothScatter(x, y, nbin = 250, raster=!TRUE, ...)
      },
      as.table=TRUE,      
      layout = c(app/2, 2, 1),
      asp = "iso",
      strip = function(..., bg) strip.default(..., bg ="#cce6ff"))
    
    id.firstpage = seq_len(app)

    ma = lapply(seq_len(nfig), function(i)
      {
        id.thispage = (i-1) * app + id.firstpage
        id.thispage = id.thispage[id.thispage <= numArrays]
        update(trobj, index.cond = list(id.thispage))
      })

    legspe = if(nchannels == 1) "where I<sub>1</sub> is the intensity of the array studied, and I<sub>2</sub> is the intensity of a \"pseudo\"-array that consists of the median across arrays." else "where I<sub>1</sub> and I<sub>2</sub> are the intensities of the two channels."
    
    legspe2 = if(class == "BeadLevelList") "The calculations are done on the summarized data obtained by using the function createBeadSummaryData from the package beadarray." else ""
    
    legend = sprintf("The figure <!-- FIG --> shows the MA plot for each array. M and A are defined as :<br>M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>)),<br>%s %s Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A. A trend in the lower range of A usually indicates that the arrays have different background intensities, this may be addressed by background correction. A trend in the upper range of A usually indicates saturation of the measurements, in mild cases, this may be addressed by non-linear normalisation (e.g. quantile normalisation).", legspe, legspe2)

    title = "MA plots"
    section = "Individual array quality"
    mamean = colMeans(abs(M), na.rm=TRUE)
    mastat = boxplot.stats(mamean)
    maout = sapply(seq_len(length(mastat$out)), function(x) which(mamean == mastat$out[x]))

    out = list("plot" = ma, "section" = section, "title" = title, "legend" = legend, "scores" = mamean, "outliers" = maout, "shape" = "rect")
    class(out) = "aqmobj.ma"
    return(out)   
  }

##Wrapper
aqm.maplot = function(dataprep, ...)
  {           
    MAplot = maplotdraw(dataprep$M, dataprep$A, dataprep$sN, dataprep$numArrays, dataprep$nchannels, dataprep$classori, ...)
    return(MAplot)         
  }
