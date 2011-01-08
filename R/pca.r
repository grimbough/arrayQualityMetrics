aqm.pca = function(x) {

  cl  = intgroupColours(x, withOpacity = FALSE)
  
  pca = prcomp(t(na.omit(x$dat)))
 
  key = cl$key
  if(!is.null(key))
    key$space = "top"
  
  pcafig = xyplot(PC2 ~ PC1 , data=as.data.frame(pca$x), pch=19, cex=1, col=cl$arrayColours,
    main = if(!is.null(cl$key)) draw.key(key = cl$key))

  legfactortip = if(length(x$intgroup)==0) " (You can indicate such a factor by colour using the 'intgroup' argument.)" else ""

  legend = sprintf("The figure <!-- FIG --> shows a scatterplot of the arrays along the first two principal components. You can use this plot to explore if the arrays cluster, and whether this is according to an intended experimental factor%s, or according to unintended causes such as &quot;batch effects&quot;. Move the mouse over the points to see the sample names.", legfactortip)
  
  new("aqmReportModule",
      plot    = pcafig,
      section = "Between array comparison",
      title   = "Principal Component Analysis",
      legend  = legend,
      shape   = list("h" = 6, "w" =6),
      svg     = if(x$usesvg)
        new("svgParameters",
              name = "pca",
              numPlotObjects = x$numArrays,
              getPlotObjNodes = getPlotPoints,
              stroke = matrix(c("1", "6", "1", "1"), nrow=2, dimnames = list(NULL, c("width", "opacity")))) else new("svgParameters")
      ) ## new
}

