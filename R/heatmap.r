#Function to perform the heatmap
hmap = function(expressionset, sN, outM, numArrays, intgroup, ...)
  {
    colourRange = rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))

    dend = as.dendrogram(hclust(outM, method = "single"))
    ord = order.dendrogram(dend)
    m = as.matrix(outM)
    colnames(m) = sN
    rownames(m) = sN
    
    if(!(missing(intgroup)||is.na(intgroup))) {
      covar  = lapply(seq(along = intgroup), function(i) pData(expressionset)[, intgroup[i]])
      lev    = lapply(seq(along = intgroup), function(i) levels(as.factor(unlist(covar[i]))))
      corres = lapply(seq(along = intgroup), function(i) matrix(0, nrow=length(lev[[i]]), ncol=2))
      colourCov = lapply(seq(along = intgroup), function(i) brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[7-i]))
      
      key = lapply(seq(along = intgroup), function(i)
        list(rect = list(col=unlist(colourCov[i])[as.factor(levels(as.factor(unlist(covar[i]))))]),
             text = list(levels(as.factor(unlist(covar[i]))))))
      key = unlist(key, recursive=FALSE)
      key$rep = FALSE
      
      foo = draw.key(key = key)

      thelegend = list(
          top = list(fun=dendrogramGrob, args=list(x=dend,side="top")),
          right=list(fun=dendrogramGrob, args=list(x=dend, side="right",
              size.add=1, add = sapply(seq(along = intgroup), function(i)
                            list(rect = list(col = "transparent",
                                             fill = unlist(colourCov[i])[as.factor(unlist(covar[i]))]))),
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
      col.regions=colourRange, main = foo, ...)
        
    leghmspe = if(is(expressionset, "BeadLevelList")) "The values used are the summarized ones obtained by using the function createBeadSummaryData from the package beadarray." else "" 

    legend = sprintf("The figure <!-- FIG --> shows a false colour heatmap of between array distances (see below for their definition). The colour scale is chosen to cover the range of distances encountered in the dataset. Arrays for which the sum of the distances to the others is much different from the others are detected as outlier arrays. The dendrogram on this plot can help to find batch effects, as well as reveal clustering of the arrays according to biological effects.<br>The distances <i>d<sub>xy</sub></i> between arrays <i>x</i> and <i>y</i> are the mean absolute difference (L<sub>1</sub>-metric) of the vector of M-values for each pair of arrays (computed on the data from all probes without filtering): <i>d<sub>xy</sub> =  mean|M<sub>xi</sub>-M<sub>yi</sub>|</i>. Here, <i>M<sub>xi</sub></i> is the M-value of the <i>i</i>-th probe on the <i>x</i>-th array. %s",  leghmspe)
## "Consider the following decomposition of <i>M<sub>xi</sub>: M<sub>xi</sub> = z<sub>i</sub> + &beta;<sub>xi</sub> + &epsilon;<sub>xi</sub></i>, where <i>z<sub>i</sub></i> is the probe effect for probe <i>i</i> (the same across all arrays), <i>&epsilon;<sub>xi</sub></i> are i.i.d. random variables with mean zero and <i>&beta;<sub>xi</sub></i> is such that for any array <i>x</i>, the majority of values <i>&beta;<sub>xi</sub></i> are negligibly small (i. e. close to zero). <i>&beta;<sub>xi</sub></i> represents differential expression effects. In this model, all values <i>d<sub>xy</sub></i> are (in expectation) the same, namely 2 times the standard deviation of <i>&epsilon;<sub>xi</sub></i>."
    
    title = "Heatmap representation of the distances between arrays"
    section = "Between array comparison"

    madsum  = rowSums(as.matrix(outM), na.rm=TRUE)
    madstat = boxplot.stats(madsum)
    madout =  which(madsum %in% madstat$out)

    out = list("plot"= hfig,
               "section" = section,
               "title" = title,
               "legend" = legend,
               "scores" = madsum,
               "outliers" = madout,
               "shape" = "square")
    class(out) = "aqmobj.heat"
    return(out) 
  }


## heatmap
aqm.heatmap = function(expressionset, dataprep, intgroup, ...)
  hmap(expressionset, dataprep$sN, dataprep$outM, dataprep$numArrays, intgroup, ...)            
