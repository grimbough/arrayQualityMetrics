## heatmap
aqm.heatmap = function(x)
{
  colourRange = rgb(seq(0, 1, l=256),
                    seq(0, 1, l=256),
                    seq(1, 0, l=256))
  
  m = as.matrix(x$outM)
  madsum  = rowSums(m, na.rm=TRUE)
  madstat = boxplot.stats(madsum)
  outliers =  which(madsum %in% madstat$out)
  
  dend = as.dendrogram(hclust(x$outM, method = "single"))
  ord = order.dendrogram(dend)

  colnames(m) = rownames(m) = paste(ifelse(seq_len(x$numArrays) %in% outliers, "* ", ""),
                                    seq_len(x$numArrays), sep="")

  if(length(x$intgroup)>0) {
    palettes = c("Set1", "Set2", "Set3", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2")
    stopifnot(all(palettes %in% rownames(brewer.pal.info)))
    palettes = rep(palettes, ceiling(length(x$intgroup) / length(palettes)))
  
    covar  = lapply(seq(along = x$intgroup), function(i) pData(x$expressionset)[[x$intgroup[i]]])
    lev    = lapply(seq(along = x$intgroup), function(i) levels(as.factor(covar[[i]])))

    colourCov = lapply(seq(along = x$intgroup), function(i)
      brewer.pal(brewer.pal.info[palettes[i], "maxcolors"], palettes[i])) 
      
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
        
  leghmspe = if(is(x$expressionset, "BeadLevelList"))
    "The values used are the summarized ones obtained by using the function createBeadSummaryData from the package beadarray." else "" 

  legend = sprintf("The figure <!-- FIG --> shows a false colour heatmap of between array distances (see below for their definition). The colour scale is chosen to cover the range of distances encountered in the dataset. Arrays for which the sum of the distances to the others is much different from the others are detected as outlier arrays (marked by &quot;*&quot;). The dendrogram on this plot can help to find batch effects, as well as reveal clustering of the arrays according to biological effects.<br>The distance <i>d<sub>xy</sub></i> between arrays <i>x</i> and <i>y</i> is the mean absolute difference (L<sub>1</sub>-distance) between the vectors of M-values of the arrays (computed on the data from all probes without filtering). In formula, <i>d<sub>xy</sub> =  mean|M<sub>xi</sub>-M<sub>yi</sub>|</i>. Here, <i>M<sub>xi</sub></i> is the M-value of the <i>i</i>-th probe on the <i>x</i>-th array. %s",  leghmspe)
## "Consider the following decomposition of <i>M<sub>xi</sub>: M<sub>xi</sub> = z<sub>i</sub> + &beta;<sub>xi</sub> + &epsilon;<sub>xi</sub></i>, where <i>z<sub>i</sub></i> is the probe effect for probe <i>i</i> (the same across all arrays), <i>&epsilon;<sub>xi</sub></i> are i.i.d. random variables with mean zero and <i>&beta;<sub>xi</sub></i> is such that for any array <i>x</i>, the majority of values <i>&beta;<sub>xi</sub></i> are negligibly small (i. e. close to zero). <i>&beta;<sub>xi</sub></i> represents differential expression effects. In this model, all values <i>d<sub>xy</sub></i> are (in expectation) the same, namely 2 times the standard deviation of <i>&epsilon;<sub>xi</sub></i>."
    
  title = "Heatmap representation of the distances between arrays"
  section = "Between array comparison"

  shape = list("h" = 6 + x$numArrays * 0.1, 
               "w" = 5 + x$numArrays * 0.1)
  
  ##-- no svg annotation: is not especially helpful, and does not scale well for large arrays
  ##  "svg" = list(annotation=annotation, getfun = function(doc) heatmapRectangles(doc, n=numArrays)))

  new("aqmReportModule",
      "plot"     = hfig,
      "section"  = section,
      "title"    = title,
      "legend"   = legend,
      "shape"    = shape)
}

## Find the numArrays x numArrays rectangles of the heatmap
##   (there is some subtlety since there are also other rectangles in the plot, from the
##    colour key, the sidebar, the legend; we avoid these by selecting 'getPlotRegionNodes')

heatmapRectangles = function(doc, n) {
  
  prn = SVGAnnotation::getPlotRegionNodes(doc)[[1]]

  ## Do some testing: are these really the n*n rectangles of the heatmap?
  ## see the man page of getBoundingBox in the SVGAnnotate package
  bb = XML::xmlSApply(prn, function(x) SVGAnnotation::getBoundingBox(x))
  
  good = ((is.matrix(bb) && (nrow(bb)==4) && (ncol(bb)==n^2) && all(colnames(bb)=="path")))
  if(good){
    ## rearrange bounding coordinates for easier comparison
    dx = array(t(bb), dim=c(n, n, 4))
    for(i in 1:2) good = good && all(diff(t(dx[,,i]))==0)
    for(i in 3:4) good = good && all(diff(  dx[,,i] )==0)
  }
  if(!good)
    stop("Error in identifying the heatmap elements in the SVG document while trying to add interactive annotation.")

  ## Now we are happy and can move on.
  return(XML::xmlChildren(prn))
}


