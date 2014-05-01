##----------------------------------------------------------------------------------
## This function is derived from boxplot.stats. The main difference is that it only
## detects outliers to the right (i.e. extremely large values).
##----------------------------------------------------------------------------------
boxplotOutliers = function (x, coef = 1.5) 
{
    stats = stats::fivenum(x, na.rm = TRUE)
    iqr = diff(stats[c(2, 4)])
    th = (stats[4] + coef * iqr)
    list(threshold = th, which = which(x > th))
}

##---------------------------------------------------------------
## Different methods for outlier detection from empiricical
## distributions
##---------------------------------------------------------------
outliers = function(exprs, method = c("KS", "sum", "upperquartile"))
{

  s = switch(method,
    KS = {
      fx = ecdf(as.vector(exprs))
      suppressWarnings(apply(exprs, 2, function(v)
        ks.test(v, y = fx, alternative="two.sided")$statistic))
    },
    sum = {
      colSums(exprs, na.rm=TRUE)
    },
    upperquartile = {
      apply(exprs, 2, quantile, na.rm=TRUE, probs=0.75)
    },
    stop(sprintf("Invalid method '%s'", method))
    )

  fo = boxplotOutliers(s)
  new("outlierDetection",
      statistic   = s,
      threshold   = fo$threshold,
      which       = fo$which)
}


aqm.outliers = function(m)
{
  values = rev(m@outliers@statistic)
  colors = rev(m@colors)
  th     = m@outliers@threshold
  n      = length(m@outliers@which)
  
  if(length(colors)==1)
    {
      colors = rep(colors, length(values))
    } else {
      stopifnot(length(colors)==length(values))  ## this should never happen
    }

  bp = function()
    {
      par(mai=c(0.6, 0.5, 0.1, 0.2))
      xlim = c(min(values, na.rm=TRUE), max(values, th, na.rm=TRUE))
      xlim = xlim + diff(xlim)*c(-1, 1)*0.04
      #
      y = seq_along(values)
      ylim = y[c(1, length(y))] + 0.5 * c(-1,1)
      #
      plot(x = values, y = y, pch = 16, cex = 1.4, col = colors,  xlab = "", ylab = "", yaxt = "n", xlim = xlim, ylim = ylim, bty = "n")
      dx = diff(par("usr")[1:2]) * 0.025
      x0 = rep(par("usr")[1], length(values))
      x1 = values
      wh = (x1 - x0) > 2*dx
      segments(x0 = x0[wh] + dx, x1 = x1[wh] - dx, y0 = y[wh], y1 = y[wh], col = colors[wh], lwd = 2)
      abline(v = th, lwd = 2)
      text(par("usr")[1], y, paste(rev(y)), adj = c(1, 0.5), xpd=NA)
     }

  mid = "exceeded the threshold and"
  legend = paste0("The figure <!-- FIG --> shows a bar chart of the ", m@outliers@description[1], 
    ", the outlier detection criterion from the previous figure. ",
    "The bars are shown in the original order of the arrays. ", 
    switch(m@outliers@description[2],
           "data-driven" = paste("Based on the distribution of the values across all arrays, a threshold of", signif(th, 3),
                                 "was determined"),
           "fixed" = paste("A threshold of", signif(th, 3), "was used"),
           stop(paste0("Invalid threshold determination method '", m@outliers@description[2], "'"))),
    ", which is indicated by the vertical line. ",
    if (n==0)
      paste("None of the arrays", mid, "was considered an outlier.") else if (n==1)
      paste("One array", mid, "was considered an outlier.") else
      paste(n, "arrays", mid, "were considered outliers."))
  
  new("aqmReportModule",
      plot    = bp,
      section = m@section,
      title   = paste("Outlier detection for", m@title),
      id      = paste("out", m@id),
      legend  = legend,
      size    = c(w = 4, h = 1 + length(values) * 0.1),
      colors  = m@colors,
      defaultdisplay = "none")
}
