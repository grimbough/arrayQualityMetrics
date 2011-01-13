##----------------------------------------------------------------------------------
## This function is derived from boxplot.stats. The main difference is that it only
## detects outliers to the right (i.e. extremely large values).
##----------------------------------------------------------------------------------
findOutliers = function (x, coef = 1.5) 
{
    nna = !is.na(x)
    stats = stats::fivenum(x, na.rm = TRUE)
    iqr = diff(stats[c(2, 4)])
    out = if (!is.na(iqr)) {
            is.na(x) | (x > (stats[4] + coef * iqr))
        }
    which(out)
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
  list(statistic = s, outliers = findOutliers(s))
}

outlierMethodExplanation = c(
  KS     = "Kolmogorov-Smirnov statistic between each array's distribution and the distribution of the pooled data",
  median = "median of these values from each array")

outlierPhrase = function(method, n)
  sprintf("Outlier detection was performed by computing the %s. For %s, this value was exceptionally large (see the manual page of <tt>findOutliers</tt> for details)%s.",
    outlierMethodExplanation[method],
    if (n==0) "none of the arrays" else if (n==1) "one array" else paste(n, "arrays"),
    if (n==0) "" else if (n==1) ", and it was marked as an outlier" else ", and they were marked as outliers")

