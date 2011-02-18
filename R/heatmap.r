## heatmap
aqm.heatmap = function(x, maxDendrogramSize = 18)
{
  colourRange = rgb(seq(0, 1, l=256),
                    seq(0, 1, l=256),
                    seq(1, 0, l=256))

  m   = dist2(x$M)
  out = outliers(m, method = "sum")
  out@description = c("sum of distances to other arrays <i>S<sub>a</sub></i>", "data-driven")
  
  dend = as.dendrogram(hclust(as.dist(m), method = "single"))
  ord = order.dendrogram(dend)

  colnames(m) = rownames(m) = paste(ifelse(seq_len(x$numArrays) %in% out@which, "* ", ""),
                                    seq_len(x$numArrays), sep="")

  ## Shall we do dendrograms?
  doDend = (ncol(m)<=maxDendrogramSize)
  thelegend = if(doDend)
    {
      list(right = list(fun=dendrogramGrob, args=list(x=dend, side="right")))
    } else {
      dummy = list()
      attr(dummy, "height") = attr(dend, "height")
      attr(dummy, "position") = numeric(0)
      list(right = list(fun=dendrogramGrob, args=list(x=dummy, ord=ord, side="right")))
    }
  
  ## Shall we draw side bars?
  ng = length(x$intgroup)
  if(ng > 0) {
    palettes = c("Set1", "Set2", "Set3", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2")
    stopifnot(all(palettes %in% rownames(brewer.pal.info)))
    palettes = rep(palettes, ceiling(ng/length(palettes))) ## make sure there are enough palettes, recycle if needed

    key = rects = vector(mode="list", length=ng)
    for(i in seq_len(ng))
      {
        fac  = as.factor(x$pData[[x$intgroup[i]]])
        pal  = brewer.pal(brewer.pal.info[palettes[i], "maxcolors"], palettes[i])
        pal  = rep(pal, ceiling(nlevels(fac)/length(pal)))[ seq_len(nlevels(fac)) ]
        cols = pal[as.integer(fac)]
        key[[i]] =  list(rect = list(col = pal),
                         text = list(levels(fac)))
        rects[[i]] = list(rect = list(col = "transparent",
                                 fill = cols,
                                 type = "rectangle"))
        if(i==1) out@colors = cols
      }
    
    key = unlist(key, recursive=FALSE)
    key$rep = FALSE
    thekey = draw.key(key = key)
    
    if(is.null(thelegend))
      {
        ## cf. dendrogramGrob
        key.layout = grid.layout(nrow = 1, ncol = ng,
            heights = unit(1, "null"), widths = unit(1, "lines"), respect = FALSE)
        thelegend = frameGrob(layout = key.layout)
        for (i in seq_len(ng)) {
            addi = rects[[i]]
            thelegend = placeGrob(thelegend,
                           rectGrob(gp = do.call(gpar, addi)), row = 1, col = i)
          }
      } else {
        thelegend$right$args = append(thelegend$right$args,
          list(size.add=1, add = rects))
      }
  } else {  
    thekey = NULL
  } ## if (ng>0) 
    
  hfig = levelplot(m[ord,ord],
    scales = list(x=list(rot=90)),
    legend = thelegend,
    colorkey = list(space ="left"),
    xlab="", ylab="",
    col.regions=colourRange, main = thekey)

  nout = length(out@which)
  
  legend = paste("The figure <!-- FIG --> shows a false colour heatmap of the distances between arrays. ",
    "The colour scale is chosen to cover the range of distances encountered in the dataset. ",
    "The presence of patterns in this plot",
    if(doDend) ", as well as the dendrogram," else "",
    " can indicate clustering of the arrays either because of intended biological or unintended experimental factors ",
    "(batch effects). ",
    "The distance <i>d<sub>ab</sub></i> between two arrays <i>a</i> and <i>b</i> is computed as the mean absolute difference ",
    "(L<sub>1</sub>-distance) between the data of the arrays (using the data from all probes without filtering). ",
    "In formula, <i>d<sub>ab</sub></i> = mean | <i>M<sub>ai</sub> - M<sub>bi</sub></i> |, ",
    "where <i>M<sub>ai</sub></i> is the value of the <i>i</i>-th probe on the <i>a</i>-th array. ",
    "Outlier detection was performed by looking for arrays for which the sum of the distances to all other arrays, ",
    "<i>S<sub>a</sub></i> = &Sigma;<sub><i>b</i></sub> <i>d<sub>ab</sub></i> was exceptionally large. ",
    if(nout>0) paste(if(nout>1) paste(nout, "such arrays were detected, and they are") else
                     "One such array was detected, and it is", "marked by an asterisk, *.") else
                        "No such arrays were detected.", sep="")
    
  new("aqmReportModule",
      plot      = hfig,
      section   = "Between array comparison",
      title     = "Distances between arrays",
      legend    = legend,
      size      = c(w = 5 + x$numArrays * 0.1, h = 6 + x$numArrays * 0.1),
      outliers  = out)
}


