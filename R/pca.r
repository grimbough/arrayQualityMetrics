aqm.pca = function(x, outliers) {

  sN = sampleNames(x$expressionset)
  annotation = namedEmptyList(x$numArrays)
  for(i in seq(along=annotation))
    annotation[[i]] = list(title = sprintf("Array %d: %s", i, sN[i]),
                           linkedids=names(annotation)[i])

  cl  = intgroupColours(x, withOpacity = FALSE)
  cex = ifelse(seq(along=sN) %in% outliers, 3, 1)
  pch = 19 
  
  pca = prcomp(t(na.omit(x$dat)))
 
  key = cl$key
  if(!is.null(key))
    key$space = "top"
  
  pcafig = xyplot(PC2 ~ PC1 , data=as.data.frame(pca$x), pch=pch, cex=cex, col=cl$arrayColours,
    main = if(!is.null(cl$key)) draw.key(key = cl$key))

  legfactortip = if(length(x$intgroup)==0) " (You can indicate such a factor by colour using the 'intgroup' argument.)" else ""
  legoutliers = if(length(outliers)>0) " Outliers -according to the same criterion as in the heatmap plot- are indicated by larger symbols." else "" 
  legend = sprintf("The figure <!-- FIG --> shows a scatterplot of the arrays along the first two principal components. You can use this plot to explore if the arrays cluster, and whether this is according to an intended experimental factor%s, or according to unintended causes such as &quot;batch effects&quot;. Move the mouse over the points to see the sample names.%s", legfactortip, legoutliers)
  
  title = "Principal Component Analysis"
  section = "Between array comparison"
   
  new("aqmTrellis",
      plot = pcafig,
      section = section,
      title = title,
      legend = legend,
      shape = list("h" = 6, "w" =6),
      svg = list(annotation=annotation, getfun=SVGAnnotation::getPlotPoints))
}

