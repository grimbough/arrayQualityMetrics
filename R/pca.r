aqm.pca = function(x, ...)
{
  pca = prcomp(t(na.omit(x$M)))

  pcafig = xyplot(PC2 ~ PC1 , data=as.data.frame(pca$x), pch=19, cex=1, col=x$arrayColors,
    main = if(!is.null(x$key)) draw.key(key = x$key), aspect = "iso")

  legend = paste("The figure <!-- FIG --> shows a scatterplot of the arrays along the first two principal components. ",
    "You can use this plot to explore if the arrays cluster, and whether this is according to an intended experimental factor",
    if(length(x$intgroup)==0) " (you can indicate such a factor by color using the 'intgroup' argument)" else "",
    ", or according to unintended causes such as batch effects. Move the mouse over the points to see the sample names.<BR>",
    "Principal component analysis is a dimension reduction and visualisation technique that is here used to project ",
    "the multivariate data vector of each array into a two-dimensional plot, such that the spatial arrangement of the ",
    "points in the plot reflects the overall data (dis)similarity between the arrays.",
    sep="")

  new("aqmReportModule",
      plot    = pcafig,
      section = "Between array comparison",
      title   = "Principal Component Analysis",
      id      = "pca",
      legend  = legend,
      size    = c(w = 1, h = 1) * 4 +  0.2 * sqrt(x$numArrays) + c( w = 0, h = 1) * length(x$key$rect$col) * 0.2,
      colors  = x$arrayColors,
      svg     = new("svgParameters",
           numPlotObjects = x$numArrays,
           getPlotObjNodes = getPlotPoints))
}

