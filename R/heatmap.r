## heatmap
aqm.heatmap = function(x)
{
  colourRange = rgb(seq(0, 1, l=256),
                    seq(0, 1, l=256),
                    seq(1, 0, l=256))

  m   = dist2(x$M)
  out = outliers(m, method = "sum")
  out@description = "sum of distances to other arrays <i>S<sub>x</sub></i>"
  
  dend = as.dendrogram(hclust(as.dist(m), method = "single"))
  ord = order.dendrogram(dend)

  colnames(m) = rownames(m) = paste(ifelse(seq_len(x$numArrays) %in% out@which, "* ", ""),
                                    seq_len(x$numArrays), sep="")

  if(length(x$intgroup)>0) {
    palettes = c("Set1", "Set2", "Set3", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2")
    stopifnot(all(palettes %in% rownames(brewer.pal.info)))
    palettes = rep(palettes, ceiling(length(x$intgroup) / length(palettes)))
  
    covar  = lapply(seq(along = x$intgroup), function(i) x$pData[[x$intgroup[i]]])
    lev    = lapply(seq(along = x$intgroup), function(i) levels(as.factor(covar[[i]])))

    colourCov = lapply(seq(along = x$intgroup), function(i)
      brewer.pal(brewer.pal.info[palettes[i], "maxcolors"], palettes[i])) 
      
    out@colors = colourCov[[1]]
    
    key = lapply(seq(along = x$intgroup), function(i) {
      fac = as.factor(covar[[i]])
      list(rect = list(col = colourCov[[i]][as.factor(levels(fac))]),
           text = list(levels(fac)))
    })
    
    key = unlist(key, recursive=FALSE)
    key$rep = FALSE
    
    foo = draw.key(key = key)
    
    thelegend = list(
      top = list(fun=dendrogramGrob, args=list(x=dend,side="top")),
      right=list(fun=dendrogramGrob, args=list(x=dend, side="right", size.add=1,
         add = sapply(seq(along = x$intgroup), function(i)
            list(rect = list(col = "transparent",
                 fill = colourCov[[i]][as.factor(covar[[i]])]))),
                                       type = "rectangle")))
  } else {
    thelegend = list(
      top  = list(fun=dendrogramGrob, args=list(x=dend, side="top")),
      right= list(fun=dendrogramGrob, args=list(x=dend, side="right")))
    foo = NULL
  }
    
  hfig = levelplot(m[ord,ord],
    scales = list(x=list(rot=90)),
    legend = thelegend,
    colorkey = list(space ="left"),
    xlab="", ylab="",
    col.regions=colourRange, main = foo)

  nout = length(out@which)
  
  legend = paste("The figure <!-- FIG --> shows a false colour heatmap of the distances between arrays.",
    "The colour scale is chosen to cover the range of distances encountered in the dataset. The dendrogram on this",
    "plot can help to find batch effects, as well as reveal clustering of the arrays according to biological effects.",
    "The distance <i>d<sub>ab</sub></i> between two arrays <i>a</i> and <i>b</i> is computed as the mean absolute difference",
    "(L<sub>1</sub>-distance) between the data of the arrays (using the data from all probes without filtering).",
    "In formula, <i>d<sub>ab</sub></i> = mean | <i>M<sub>ai</sub> - M<sub>bi</sub></i> |,",
    "where <i>M<sub>ai</sub></i> is the value of the <i>i</i>-th probe on the <i>a</i>-th array.",
    "Outlier detection was performed by looking for arrays for which the sum of the distances to all other arrays, ",
    "<i>S<sub>a</sub></i> = &Sigma;<sub><i>b</i></sub> <i>d<sub>ab</sub></i> was exceptionally large.",
    if(nout>0) paste(if(nout>1) paste(nout, "such arrays were detected, and they are") else
                     "One such array was detected, and it is", "marked by an asterisk, *.") else
                        "No such arrays were detected.") 
    
  new("aqmReportModule",
      plot      = hfig,
      section   = "Between array comparison",
      title     = "Distances between arrays",
      legend    = legend,
      size      = c(w = 5 + x$numArrays * 0.1, h = 6 + x$numArrays * 0.1),
      outliers  = out)
}


