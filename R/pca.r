aqm.pca = function(x) {

  cl  = intgroupColours(x, withOpacity = FALSE)
  
  pca = prcomp(t(na.omit(x$M)))
 
  key = cl$key
  if(!is.null(key))
    key$space = "top"
  
  pcafig = xyplot(PC2 ~ PC1 , data=as.data.frame(pca$x), pch=19, cex=1, col=cl$arrayColours,
    main = if(!is.null(cl$key)) draw.key(key = cl$key))

  legend = paste("The figure <!-- FIG --> shows a scatterplot of the arrays along the first two principal components. ",
    "You can use this plot to explore if the arrays cluster, and whether this is according to an intended experimental factor",
    if(length(x$intgroup)==0) " (you can indicate such a factor by colour using the 'intgroup' argument)" else "",
    ", or according to unintended causes such as batch effects. Move the mouse over the points to see the sample names.<BR>",
    "Principal component analysis is a dimension reduction and visualisation technique that is here used to project ",
    "the multivariate data vector of each array into a two-dimensional plot, such that the spatial arrangement of the ",
    "points in the plot reflects the overall data (dis)similarity between the arrays.",
    sep="")
  
  new("aqmReportModule",
      plot    = pcafig,
      section = "Between array comparison",
      title   = "Principal Component Analysis",
      legend  = legend,
      size    = c(w = 6, h =6),
      svg     = new("svgParameters",
           name = "pca",
           numPlotObjects = x$numArrays,
           getPlotObjNodes = getPlotPoints,
           stroke = matrix(c("1", "6", "1", "1"), nrow=2, dimnames = list(NULL, c("width", "opacity")))))
}

