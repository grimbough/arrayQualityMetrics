##---------------------------------------------------------------
## This function is derived from boxplot.stats. The main difference is that it only
## detects outliers to the right (i.e. extremely large values).
##---------------------------------------------------------------
outliers = function (x, coef = 1.5) 
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
## Simple distribution-based outlier detection using KS-statistic
##  (not: the p-value)
##---------------------------------------------------------------
ksOutliers = function(exprs, subsamp = 300)
{
  if (nrow(exprs)>subsamp)
    exprs = exprs[sample(nrow(exprs), subsamp), ]

  fx = ecdf(as.vector(exprs))
  
  s = suppressWarnings(apply(exprs, 2, function(v)
    ks.test(v, y = fx, alternative="two.sided")$statistic))

  list(statistic = s, outliers = outliers(s))
}

