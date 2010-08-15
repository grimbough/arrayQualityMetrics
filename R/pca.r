aqm.pca = function(expressionset, dataprep, intgroup, ...)
{

  sN = sampleNames(expressionset)
  annotation = vector(mode="list", length = length(sN))
  names(annotation) = sprintf("line%d", seq(along=annotation))
  for(i in seq(along=annotation))
    annotation[[i]] = list(title = sprintf("Array %d: %s", i, sN[i]),
                           linkedids=names(annotation)[i])

  cl = intgroupColours(intgroup, expressionset)

  pca = prcomp(t(na.omit(dataprep$dat)))
 
  key =cl$key
  key$space = "top"
  
  pcafig = xyplot(PC2 ~ PC1 , as.data.frame(pca$x), pch=19, col=cl$arrayColours, key = key)

  legend = "The figure <!-- FIG --> shows a scatterplot of the arrays along the first two principal components. We expect the arrays to cluster according to a relevant experimental factor. "
  title = "Principal Component Analysis"
  section = "Between array comparison"
   
  out = list("plot" = pcafig,
             "section" = section,
             "title" = title,
             "legend" = legend,
             "shape" = "square",
             "svg" = list(annotation=annotation, getfun=SVGAnnotation::getPlotPoints))
  
  class(out) = "aqmobj.pca"
  return(out) 
}

