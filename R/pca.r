aqm.pca = function(expressionset, dataprep, intgroup, ...)
{

  sN = sampleNames(expressionset)
  annotation = namedEmptyList(length(sN))
  for(i in seq(along=annotation))
    annotation[[i]] = list(title = sprintf("Array %d: %s", i, sN[i]),
                           linkedids=names(annotation)[i])

  cl = intgroupColours(intgroup, expressionset)

  pca = prcomp(t(na.omit(dataprep$dat)))
 
  key =cl$key
  key$space = "top"
  
  pcafig = xyplot(PC2 ~ PC1 , as.data.frame(pca$x), pch=19, col=cl$arrayColours, key = key)

  legend = "The figure <!-- FIG --> shows a scatterplot of the arrays along the first two principal components. You can use this plot to explore if the arrays cluster, and whether this is according to an intended experimental factor or according to unintended causes such as &quot;batch effects&quot;. Move the mouse over the points to see the sample names."
  title = "Principal Component Analysis"
  section = "Between array comparison"
   
  out = list("plot" = pcafig,
             "section" = section,
             "title" = title,
             "legend" = legend,
             "shape" = list("h" = 6, "w" =6),
             "svg" = list(annotation=annotation, getfun=SVGAnnotation::getPlotPoints))
  
  class(out) = "aqmobj.pca"
  return(out) 
}

