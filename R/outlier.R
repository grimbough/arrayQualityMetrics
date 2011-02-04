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
outliers = function(exprs, method = c("KS", "mean", "median"))
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
    median = {
      apply(exprs, 2, median, na.rm=TRUE)
    },
    stop(sprintf("Invalid method '%s'", method))
    )

  fo = findOutliers(s)
  new("outlierDetection",
      statistic = s,
      threshold = fo$threshold,
      which     = fo$which)
}

outlierMethodExplanation = c(
  KS     = "Kolmogorov-Smirnov statistic between each array's distribution and the distribution of the pooled data.",
  median = "median of these values from each array.")

outlierPhrase = function(method=FALSE, n)
  {
    rv = if(!identical(method, FALSE))
      {
        stopifnot(method%in%names(outlierMethodExplanation))
        paste("Outlier detection was performed by computing the", outlierMethodExplanation[method])
      } else {
        ""
      }

    rv = paste(rv, sprintf("For %s, this value was exceptionally large (see the manual page of <tt>findOutliers</tt> for details)%s.",
      if (n==0) "none of the arrays" else if (n==1) "one array" else paste(n, "arrays"),
      if (n==0) "" else if (n==1) ", and it was marked as an outlier" else ", and they were marked as outliers"))

    return(rv)
  }

