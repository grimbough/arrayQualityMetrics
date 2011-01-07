##---------------------------------------------------------------
## Simple distribution-based outlier detection using KS-statistic
##  (not: the p-value)
##---------------------------------------------------------------
ksOutliers = function(x, subsamp = 300, theta = 2)
{
  if (nrow(x)>subsamp)
    x = x[sample(nrow(x), subsamp), ]

  fx = ecdf(as.vector(x))
  
  s = suppressWarnings(apply(x, 2, function(v)
    ks.test(v, y = fx, alternative="two.sided")$statistic))

  list(statistic = s,
       outliers = which( (s-mean(s)) / sd(s) > theta ))
}


##----------------------------------------
## aqm.boxplot
##----------------------------------------
aqm.boxplot = function(x, subsample = 10000) {
  
  ks = ksOutliers(x$dat)
  
  if (nrow(x$dat)>subsample) {
    ss = sample(nrow(x$dat), subsample)
    nr = length(ss)
  } else {
    ss = TRUE
    nr = nrow(x$dat)
  }
  
  sample_id = rep( seq_len(x$numArrays), each = nr )
  
  if(x$nchannels == 2)  {
    values    = with(x, c(rc[ss,], gc[ss,], dat[ss,]))
    sample_id = rep(sample_id, times = 3)
    panels    = factor(rep(1:3, each = nr * x$numArrays),
                   levels = 1:3,
                   labels = c("a. Red Channel", "b. Green Channel", "c. Log2(Ratio)"))
    formula = sample_id ~ values | panels
    lay = c(3,1)
    legspe = "Three panels are shown: left, red channel; middle, green channel; right, log<sub>2</sub>(ratio). Outlier detection  was performed on the distribution of log<sub>2</sub>(ratio). "
   } else {
    values  = as.numeric(x$dat[ss, ])
    formula = sample_id ~ values
    lay = c(1,1)
    legspe = ""
  }
  xAsterisk = quantile(x$dat[ss,], probs=0.01, na.rm=TRUE)
  
  cl = intgroupColours(x)

  box = bwplot(formula, groups = sample_id, layout = lay, as.table = TRUE,
        strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
        horizontal = TRUE,
        main = if(!is.null(cl$key)) draw.key(key = cl$key), 
        pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
        xlab = "", ylab = "Array",
        fill = cl$arrayColours, panel = panel.superpose, 
        scales = list(x=list(relation="free"), y=list(axs="i")),
        ylim = c(x$numArrays+0.7,0.3),
        prepanel = function(x, y) {
          list(xlim = quantile(x, probs = c(0.01, 0.99), na.rm=TRUE))
        },
        panel.groups = function(x, y, ...) {
          panel.bwplot(x, y, ...)
          if(lattice:::packet.number()==lay[1]) {
            whArray = list(...)$group.value
            if (whArray %in% ks$outliers)
              ltext(xAsterisk, whArray, "*", font=2, cex=2, adj=c(0.5, 0.75))
          }
        })

  shape = list("h" = 2.5 + x$numArrays * 0.1 +  1/x$numArrays, 
               "w" = 3+3*lay[1])
  
  outliertext = if(length(ks$outliers)>0) "Outliers -according to the Kolmogorov-Smirnov statistic between each array's distribution and the distribution of the pooled data- are marked by an asterisk (*). " else ""
  
  legend = sprintf("The figure <!-- FIG --> presents boxplots representing summaries of the signal intensity distributions of the arrays. %sEach box corresponds to one array. Typically, one expects the boxes to have similar positions and widths. If the distribution of an array is very different from the others, this may indicate an experimental problem. %s", legspe, outliertext)

  title = "Boxplots"
  section = "Array intensity distributions"

  new("aqmReportModule",
      plot = box,
      section = section,
      title = title,
      legend = legend,
      outliers = ks$outliers,
      shape = shape)
}
