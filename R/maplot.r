aqm.maplot = function(x, subsample=20000, Dthresh=0.15) {

  if(x$nchannels==1)
    {
      stopifnot(identical(x$M, x$A))
      medArray = rowMedians(x$M, na.rm=TRUE)
      M =  x$M - medArray
      A = (x$M + medArray)/2
    } else {
      M = x$M
      A = x$A
      stopifnot(identical(dim(M), dim(A)))
    }

  ## For each array j, compute the D statistic from Hoeffding's test for independence 
  ## and sort / detect outlier arrays by the value of the test statistic. We do subsampling
  ## to save time
  if(nrow(M)>subsample)
    {
      sel = sample(nrow(M), subsample)
      sM = M[sel, ]
      sA = A[sel, ]
    } else {
      sM = M
      sA = A
    }
  stat = sapply(seq_len(x$numArrays), function(j)
    {
      hoeffd(sA[,j], sM[,j])$D[1,2]
    })
  maout = which(stat > Dthresh)
    
  ## Plot maximally 8 scatterplots
  if(x$numArrays<=8)
    {
      whj = seq_len(x$numArrays)
      lay = c(ceiling(x$numArrays/2), 2)
      legOrder = ""
    } else {
      whj = order(stat, decreasing=TRUE)[c(1:4, x$numArrays+c(-3:0))]
      lay = c(4, 2)
      legOrder = "Shown are the 4 arrays with the highest value of <i>D</i> (top row), and the 4 arrays with the lowest value (bottom row). "
    }
    
  xlim = quantile(A, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
  ylim = quantile(M, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
  panelNames = sprintf("array %d (D=%4.2f)", whj, stat[whj]) 

  i = seq(along=whj)
  df = data.frame(
    i = factor(i, levels = i),
    px = i,
    py = i)
  
  ma = xyplot(py ~ px | i,
    data = df,
    xlim = xlim,
    ylim = ylim,
    xlab = "A",
    ylab = "M",
    panel = function(x, y, ...) panel.smoothScatter(x=A[, whj[x]], y=M[, whj[y]], raster=TRUE, nbin=250, ...),
    as.table = TRUE,      
    layout = lay,
    asp = "iso",
    strip = function(..., bg, factor.levels) strip.default(..., bg ="#cce6ff", factor.levels = panelNames))
  
  legRef = if(x$nchannels == 1)
    "where I<sub>1</sub> is the intensity of the array studied, and I<sub>2</sub> is the intensity of a \"pseudo\"-array that consists of the median across arrays." else
    "where I<sub>1</sub> and I<sub>2</sub> are the intensities of the two channels."
    
  legend = sprintf("The figure <!-- FIG --> shows the MA plot for each array. M and A are defined as :<br>M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>)),<br>%s Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, and there should be no trend in M as a function of A. If there is a trend in the lower range of A, this often indicates that the arrays have different background intensities; this may be addressed by background correction. A trend in the upper range of A can indicate saturation of the measurements; in mild cases, this may be addressed by non-linear normalisation (e.g. quantile normalisation).<br>Outlier detection was performed by computing Hoeffding's <i>D</i>-statistic on the joint distribution of A and M for each array. %sThe value of <i>D</i> is shown in the panel headings. %s had <i>D</i>&gt;%g and %s marked as outliers. For more information on Hoeffing's <i>D</i> statistic, please see the manual page of the function <tt>hoeffd</tt> in the <tt>Hmisc</tt> package.", legRef, legOrder, if(length(maout)==1) "One array" else paste(length(maout), "arrays"), Dthresh, if(length(maout)==1) "was" else "were")
  
  new("aqmReportModule",
      plot = ma, 
      section = "Individual array quality",
      title = "MA plots",
      legend = legend,
      outliers = maout,
      shape = list(h=6, w=10))
}

