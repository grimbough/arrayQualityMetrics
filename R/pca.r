aqm.pca = function(expressionset, dataprep, intgroup = "Covariate", ...)
{
    covar = pData(expressionset)[colnames(pData(expressionset))==intgroup[1]][,1]
    colourCovd = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[6])

    cols = unique(colourCovd[as.factor(covar)])

    pca = prcomp(t(dataprep$dat))
 
    key = list(points = list(pch = 19, col=cols), text = list(levels(as.factor(unlist(covar)))))
    key$rep = FALSE
    key$space = "right"
    pcafig = xyplot(PC2 ~ PC1 , as.data.frame(pca$x), pch=19, col=cols, key = key)

    legend = "Figure <!-- FIG --> represents the Principal Component Analysis of the arrays. The colours correspond to the group of interest given."
    title = "Principal Component Analysis of the arrays"
    section = "Between array comparison"
   
    out = list("plot" = pcafig, "section" = section,  "title" = title, "legend" = legend, "shape" = "square")
    class(out) = "aqmobj.pca"
    return(out) 
}

