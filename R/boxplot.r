ksOutliers = function(x, subsamp = 300, theta = 2){
  if (nrow(x)>subsamp)
    x = x[sample(nrow(x), subsamp), ] 
  s = apply(x, 2, function(v)
    suppressWarnings(ks.test(v, x, alternative="two.sided")$statistic))
  which( (s-mean(s)) / sd(s) > theta )
}


aqm.boxplot = function(expressionset, dataprep, intgroup, outliers, subsample = 10000, ...) {

  if (nrow(dataprep$dat)>subsample) {
    ss = sample(nrow(dataprep$dat), subsample)
    nr = length(ss)
  } else {
    ss = TRUE
    nr = nrow(dataprep$dat)
  }
  
  sample_id = rep( seq_len(dataprep$numArrays), each = nr )
  
  if(dataprep$nchannels == 2)  {
    values    = with(dataprep, c(rc[ss,], gc[ss,], dat[ss,]))
    sample_id = rep(sample_id, times = 3)
    panels    = factor(rep(1:3, each = nr * dataprep$numArrays),
                   levels = 1:3,
                   labels = c("a. Red Channel", "b. Green Channel", "c. Log(Ratio)"))
    formula = sample_id ~ values | panels
    lay = c(3,1)
  } else {
    values    = with(dataprep, dat[ss, ])
    formula = sample_id ~ values
    lay = c(1,1)
  }

  if (!missing(intgroup)) {
    groups  = as.factor(pData(expressionset)[, intgroup])
    igroups = as.integer(groups)
    colours = brewer.pal(9, "Set1")
    if(nlevels(groups) > length(colours)) {
      warning(sprintf("'intgroup' has %d levels, but only the first 9 are used for colouring.", nlevels(groups)))
      igroups[igroups > length(colours)] = length(colours)+1
      colours = c(colours, "#101010")
    } else {
      colours = colours[1:nlevels(groups)]
    }
    colsb = colours[igroups]
    foo = draw.key(key = list(
                     rect = list(col = colours),
                     text = list(levels(groups)),
                     rep = FALSE))
    
  } else {
    colsb = rep("#1F78B4", dataprep$numArrays)
    foo = NULL
  }

  if(!missing(outliers))
    colsb[outliers] = "grey"

  box = bwplot(formula, groups = sample_id, layout=lay, as.table=TRUE,
        strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
        horizontal = TRUE, 
        pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
        xlab = "", ylab = "Array",
        fill = colsb, panel = panel.superpose, panel.groups = panel.bwplot, main=foo, 
        scales = "free", prepanel =
           function(x, y) {
             list(xlim=quantile(x, probs = c(0.01, 0.99)),
                  ylim=c(dataprep$numArrays+0.7, 0.3))
           },    ...)


  shape = list("h" = if(dataprep$numArrays > 50) (dataprep$numArrays/8) else 10, "w" = 10) 
  
  legspe = if(dataprep$nchannels == 2) "Left panel: red channel, middle panel: green channel, right panel: log<sub>2</sub>(ratio)." else ""
  
  legend = sprintf("The figure <!-- FIG --> presents boxplots of the log<sub>2</sub>(Intensities). Each box corresponds to one array. %s It gives a simple summary of the distribution of probe intensities across all arrays. Typically, one expects the boxes to have similar size (IQR) and y position (median). If the distribution of an individual array is very different from the others, this may indicate an experimental problem. After normalisation, the distributions should be similar.", legspe)

  title = "Boxplots"
  section = "Array intensity distributions"

  out = list("plot" = box,
             "section" = section,
             "title" = title,
             "legend" = legend,
             "shape" = shape)
  
  class(out) = "aqmobj.box"
  
  return(out)   
}
