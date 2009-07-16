aqm.pca = function(expressionset, dataprep, intgroup = "Covariate", ...)
{
    covar = pData(expressionset)[colnames(pData(expressionset))==intgroup[1]][,1]
    colourCovd = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[6])

    cols = colourCovd[as.factor(covar)]

    pca = prcomp(t(dataprep$dat))
 
    key = list(points = list(pch = 19, col=unique(cols)), text = list(levels(as.factor(unlist(covar)))))
    key$rep = FALSE
    key$space = "top"
    pcafig = xyplot(PC2 ~ PC1 , as.data.frame(pca$x), pch=19, col=cols, key = key)

    legend = "The figure <!-- FIG --> represents a biplot for the first two principal components from the dataset. The colours correspond to the group of interest given. We expect the arrays to cluster accordingly to a relevant experimental factor. The principal components transformation of a data matrix re-expresses the features using linear combination of the original variables. The first principal component is the linear combination chosen to possess maximal variance, the second is the linear combination orthogonal to the first possessing maximal variance among all orthogonal combination."
    title = "Principal Component Analysis"
    section = "Between array comparison"
   
    out = list("plot" = pcafig, "section" = section,  "title" = title, "legend" = legend, "shape" = "square")
    class(out) = "aqmobj.pca"
    return(out) 
}

