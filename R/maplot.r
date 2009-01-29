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
        panel.smoothScatter(x, y ,...)
      },
      as.table=TRUE,      
      layout = c(app/2, 2, 1),
      asp = "iso",
      strip = function(..., bg) strip.default(..., bg ="#cce6ff"))
    
    id.firstpage = seq_len(app)

    ma = lapply(seq_len(nfig), function(i) {id.thispage = (i-1) * app + id.firstpage; id.thispage = id.thispage[id.thispage <= numArrays]; update(trobj, index.cond = list(id.thispage))})

    if(nchannels == 1)
      legspe = "where I<sub>1</sub> is the intensity of the array studied and I<sub>2</sub> is the intensity of a \"pseudo\"-array, which have the median values of all the arrays."
    if(nchannels == 2)
      legspe = "where I<sub>1</sub> and I<sub>2</sub> are the vectors of intensities of the two channels."
    
    legspe2 = if(class == "BeadLevelList") "The calculations are done on the summarized data obtained by using the function createBeadSummaryData from the package beadarray." else ""
    
    legend = sprintf("The <b>figure <!-- FIG --></b> represents the MA plot for each array. M and A are defined as :<br>M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>)),<br>%s %s Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A. Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).", legspe, legspe2)

    title = "MA plots"
    type = "Individual array quality"
    mamean = colMeans(abs(M), na.rm=TRUE)
    mastat = boxplot.stats(mamean)
    maout = sapply(seq_len(length(mastat$out)), function(x) which(mamean == mastat$out[x]))

    out = list("plot" = ma, "type" = type, "title" = title, "legend" = legend, "scores" = mamean, "outliers" = maout, "shape" = "rect")
    class(out) = "aqmobj.ma"
    return(out)   
  }

##Wrapper
aqm.maplot = function(dataprep, ...)
  {           
    MAplot = maplotdraw(dataprep$M, dataprep$A, dataprep$sN, dataprep$numArrays, dataprep$nchannels, dataprep$classori, ...)
    return(MAplot)         
  }
