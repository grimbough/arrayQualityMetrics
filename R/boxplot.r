##----------------------------------------
## aqm.boxplot
##----------------------------------------
aqm.boxplot = function(x, subsample=20000, outlierMethod = "KS") {
  
  if (nrow(x$M)>subsample) {
    ss  = sample(nrow(x$M), subsample)
    Mss = x$M[ss,,drop=FALSE]
  } else {
    ss  = TRUE
    Mss = x$M
  }
  out = outliers(Mss, method = outlierMethod)
  out@description = c("Kolmogorov-Smirnov statistic <i>K<sub>a</sub></i>", "data-driven")
  out@colors = x$arrayColors
    
  sample_id = rep( seq_len(x$numArrays), each = nrow(Mss) )
  
  if(x$nchannels == 2)  {
    values    = with(x, c(R[ss,], G[ss,], Mss))
    sample_id = rep(sample_id, times = 3)
    panels    = factor(rep(1:3, each = nrow(Mss) * x$numArrays),
                   levels = 1:3,
                   labels = c("a. Red Channel", "b. Green Channel", "c. Log2(Ratio)"))
    formula = sample_id ~ values | panels
    lay = c(3,1)
    legPanels = c("Three panels are shown: left, red channel; middle, green channel; right, log<sub>2</sub>(ratio). ",
                  "on the distribution of Log<sub>2</sub>(Ratio). ")
   } else {
    values  = as.numeric(Mss)
    formula = sample_id ~ values
    lay = c(1,1)
    legPanels = c("", "")
  }
  xAsterisk = quantile(Mss, probs=0.01, na.rm=TRUE)
  
  box = bwplot(formula, groups = sample_id, layout = lay, as.table = TRUE,
        strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
        horizontal = TRUE,
        main = if(!is.null(x$key)) draw.key(key = x$key), 
        pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
        xlab = "", ylab = "Array",
        fill = x$arrayColors, panel = panel.superpose, 
        scales = list(x=list(relation="free"), y=list(axs="i")),
        ylim = c(x$numArrays+0.7,0.3),
        prepanel = function(x, y) {
          list(xlim = quantile(x, probs = c(0.01, 0.99), na.rm=TRUE))
        },
        panel.groups = function(x, y, ...) {
          panel.bwplot(x, y, ...)
          if(lattice:::packet.number()==lay[1]) {
            whArray = list(...)$group.value
            if (whArray %in% out@which)
              ltext(xAsterisk, whArray, "*", font=2, cex=2, adj=c(0.5, 0.75))
          }
        })
  
  legend = paste("The figure <!-- FIG --> shows boxplots representing summaries of the signal intensity distributions of the arrays. ",
    legPanels[1],
    "Each box corresponds to one array. Typically, one expects the boxes to have similar positions and widths. If the ",
    "distribution of an array is very different from the others, this may indicate an experimental problem. ", 
    "Outlier detection was performed ",
    legPanels[2],
    "by computing the Kolmogorov-Smirnov statistic <i>K<sub>a</sub></i> between each array's distribution ",
    "and the distribution of the pooled data.",
    sep="")

  new("aqmReportModule",
      plot = box,
      section = "Array intensity distributions",
      title = "Boxplots",
      id = "box",
      legend = legend,
      outliers = out,
      colors  = x$arrayColors,
      size = c(w = 3+3*lay[1], h = 2.5 + x$numArrays * 0.1 +  1/x$numArrays + length(x$key$rect$col) * 0.2))
}
