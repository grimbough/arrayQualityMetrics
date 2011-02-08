##----------------------------------------------------------------------------------
## This function is derived from boxplot.stats. The main difference is that it only
## detects outliers to the right (i.e. extremely large values).
##----------------------------------------------------------------------------------
findOutliers = function (x, coef = 1.5) 
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
outliers = function(exprs, method = c("KS", "sum", "median"))
{

  s = switch(method,
    KS = {
      description = "Kolmogorov-Smirnov statistic between each array's distribution and the distribution of the pooled data"
      fx = ecdf(as.vector(exprs))
      suppressWarnings(apply(exprs, 2, function(v)
        ks.test(v, y = fx, alternative="two.sided")$statistic))
    },
    sum = {
      description = "sum of the values from each array"
      colSums(exprs, na.rm=TRUE)
    },
    median = {
      description = "median of the values from each array"
      apply(exprs, 2, median, na.rm=TRUE)
    },
    stop(sprintf("Invalid method '%s'", method))
    )

  fo = findOutliers(s)
  new("outlierDetection",
      statistic   = s,
      threshold   = fo$threshold,
      which       = fo$which,
      description = description)
}


aqm.outliers = function(m)
{
  values = rev(m@outliers@statistic)
  colors = rev(m@outliers@colors)
  th     = m@outliers@threshold
  n      = length(m@outliers@which)
  
  xlim = c(0, max(values, th, na.rm=TRUE))
  bp = function()
    {
      par(mai=c(0.6, 0.5, 0.1, 0.2))
      b = barplot(values, col = colors, xaxs = "r", names.arg = "",
              xlab = "", ylab = "", horiz = TRUE, xlim = xlim)
      abline(v = th, lwd = 2)
      text(par("usr")[1], b, paste(rev(seq(along=values))), adj = c(1, 0.5), xpd=NA, cex=0.66) 
    }

  mid = "exceeded the threshold and"
  legend = paste("The figure <!-- FIG --> shows a bar chart of the ", m@outliers@description, 
    ", the outlier detection criterion from the previous figure. ",
    "The bars are shown in the original order of the arrays. ", 
    "Based on the distribution of the values across all arrays, a threshold of ", signif(th, 3),
    " was determined, which is indicated by the vertical line. ",
    if (n==0)
      paste("None of the arrays", mid, "was considered an outlier.") else if (n==1)
      paste("One array", mid, "was considered an outlier.") else
      paste(n, "arrays", mid, "were considered outliers."),
    sep="")
  
  new("aqmReportModule",
      plot    = bp,
      section = m@section,
      legend  = legend,
      size    = c(w = 4, h = 1 + length(values) * 0.1))
}
