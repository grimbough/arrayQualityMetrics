## heatmap
aqm.heatmap = function(x)
{
  colorRange = rgb(seq(0, 1, l=256),
                    seq(0, 1, l=256),
                    seq(1, 0, l=256))

  m   = dist2(x$M)
  out = outliers(m, method = "sum")
  out@description = c("sum of distances to other arrays <i>S<sub>a</sub></i>", "data-driven")
  
  dend = as.dendrogram(hclust(as.dist(m), method = "single"))
  ord = order.dendrogram(dend)

  colnames(m) = rownames(m) = paste(ifelse(seq_len(x$numArrays) %in% out@which, "* ", ""),
                                    seq_len(x$numArrays), sep="")

  ## Shall we draw a dendrogram?
  if (ncol(m) <= arrayQualityMetricsGlobalParameters$maxNumberOfArraysForDrawingDendrogram)
    {
      thelegend = list(right = list(fun=dendrogramGrob, args=list(x=dend, side="right")))
    } else {
      thelegend = NULL
    }

  ## Shall we draw side bars?
  maxNrColors = 0
  ng = length(x$intgroup)
  if(ng > 0) {
    palettes = c("Set1", "Set2", "Set3", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2")
    stopifnot(all(palettes %in% rownames(brewer.pal.info)))
    palettes = rep(palettes, ceiling(ng/length(palettes))) ## make sure there are enough palettes, recycle if needed

    key = rects = vector(mode="list", length=ng)
    names(rects) = rep("rect", ng)
    
     for(i in seq_len(ng))
      {
        colors = brewer.pal(brewer.pal.info[palettes[i], "maxcolors"], palettes[i])
        fac  = as.factor(x$pData[[x$intgroup[i]]])
        fac  = maximumLevels(fac, n = length(colors)) ## make sure that factor has at most n levels
        colors = colors[seq_len(nlevels(fac))]
        ac = colors[as.integer(fac)]

        key[[i]] =  list(rect = list(col = colors),
                         text = list(levels(fac)))
        rects[[i]] = list(col = "transparent",
                          fill = ac[ord],
                          type = "rectangle")
        if(length(colors)>maxNrColors)
          maxNrColors = length(colors)
      }

      
    key = unlist(key, recursive=FALSE)
    key$rep = FALSE
    thekey = draw.key(key = key) 

    if(!is.null(thelegend))
      {
        thelegend$right$args = append(thelegend$right$args,
          list(size.add=1, add = rects))
      } else {
        ## adapted from 'latticeExtra:dendrogramGrob'
        lay= grid.layout(nrow = 1, ncol = ng,
            heights = unit(1, "null"),
            widths = unit(rep(1, length = ng), rep("lines", ng)), respect = FALSE)
        g = frameGrob(layout = lay)
        dy = 1/x$numArrays
        y  = seq_len(x$numArrays)*dy 
        for (i in seq_len(ng))
          {
            g = placeGrob(g,
              rectGrob(y = y, height = dy, vjust =1,
                       gp = do.call(gpar, rects[[i]])), row = 1, col = i)
          }
        idem = function(x) x
        thelegend = list(right=list(fun=idem, args=list(x=g)))
      }
    
  } else {  
    thekey = NULL
  } ## if (ng>0) 

  hfig = levelplot(m[ord,ord],
    scales = list(x=list(rot=90)),
    legend = thelegend,
    colorkey = list(space ="left"),
    xlab = "",
    ylab = "",
    col.regions = colorRange,
    main = thekey)

  nout = length(out@which)
  
  legend = paste("The figure <!-- FIG --> shows a false color heatmap of the distances between arrays. ",
    "The color scale is chosen to cover the range of distances encountered in the dataset. ",
    "Patterns in this plot can indicate clustering of the arrays either because of intended biological or ",
    "unintended experimental factors (batch effects). ",
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
      id        = "hm",
      legend    = legend,
      size      = c(w = 5 + x$numArrays * 0.075, h = 3 + x$numArrays * 0.075 + maxNrColors * 0.2),
      colors    = x$arrayColors,
      outliers  = out)
}


