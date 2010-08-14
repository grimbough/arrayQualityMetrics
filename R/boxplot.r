ksOutliers = function(x, subsamp = 300, theta = 2){
  if (nrow(x)>subsamp)
    x = x[sample(nrow(x), subsamp), ] 
  s = apply(x, 2, function(v)
    suppressWarnings(ks.test(v, x, alternative="two.sided")$statistic))
  which( (s-mean(s)) / sd(s) > theta )
}

##----------------------------------------
## This function returns a list with:
##   arrayColors:  color code for each array
##   key: a key explaining the mapping of factor values to colours
##----------------------------------------
intgroupColours = function(intgroup, expressionset){
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
    list(arrayColours = colours[igroups],
         key = draw.key(key = list(
                          rect = list(col = colours),
                          text = list(levels(groups)),
                          rep = FALSE)))
    
  } else {
    list(arrayColours = rep("#1F78B4", dataprep$numArrays),
         key = NULL)
  }
}

##----------------------------------------
## aqm.boxplot
##----------------------------------------
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
                   labels = c("a. Red Channel", "b. Green Channel", "c. Log2(Ratio)"))
    formula = sample_id ~ values | panels
    lay = c(3,1)
   } else {
    values  = dataprep$dat[ss, ]
    formula = sample_id ~ values
    lay = c(1,1)
  }
  xAsterisk = quantile(dataprep$dat[ss,], probs = 0.01)
  
  cl = intgroupColours(intgroup, expressionset)

  if(missing(outliers))
    outliers=c()
  
  box = bwplot(formula, groups = sample_id, layout = lay, as.table = TRUE,
        strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
        horizontal = TRUE, main = cl$key, 
        pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
        xlab = "", ylab = "Array",
        fill = cl$arrayColours, panel = panel.superpose, 
        scales = list(x=list(relation="free"), y=list(axs="i")),
        ylim = c(dataprep$numArrays+0.7,0.3),
        prepanel = function(x, y) {
          list(xlim = quantile(x, probs = c(0.01, 0.99)))
        },
        panel.groups = function(x, y, ...) {
          panel.bwplot(x, y, ...)
          if(packet.number()==lay[1]) {
            whArray = list(...)$group.value
            if (whArray %in% outliers)
              ltext(xAsterisk, whArray, "*", font=2, cex=3, adj=c(0.5,0.75))
          }
        },
    ...)


  shape = list("h" = 2.5 + dataprep$numArrays * 0.1 +  1/dataprep$numArrays, 
               "w" = 2.5*lay[1] + 2.5)
  
  legspe = if(dataprep$nchannels == 2) "Left panel: red channel, middle panel: green channel, right panel: log<sub>2</sub>(ratio)." else ""
  
  legend = sprintf("The figure <!-- FIG --> presents boxplots of the data. It gives a simple summary of the signal intensity distributions across all arrays. %s Each box corresponds to one array. Typically, one expects the boxes to have similar positions and widths. If the distribution of an individual array is very different from the others, this may indicate an experimental problem. After normalisation, the distributions should be more similar.", legspe)

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
