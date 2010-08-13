setClassUnion("aqmInputObj", c("ExpressionSet", "AffyBatch", "NChannelSet", "BeadLevelList"))

## set the method for arrayQualityMetrics
setGeneric("arrayQualityMetrics",
           function(expressionset,
                    outdir = getwd(),
                    force = FALSE,
                    do.logtransform = FALSE,
                    intgroup = "Covariate",
                    grouprep = FALSE,
		    spatial = TRUE,
                    sN = NULL)
           standardGeneric("arrayQualityMetrics"))

## aqmInputObj
setMethod("arrayQualityMetrics", signature(expressionset = "aqmInputObj"),
  function(expressionset, outdir, force, do.logtransform, intgroup, grouprep, spatial, sN) {

    ## Argument checking: 
    ## (Done here, once and for all, rather than where the arguments are actually consumed - 
    ## in the hope that this simplifies the code and the user interface)
    for(v in c("outdir", "intgroup"))
      if (!(is.character(get(v)) && (length(get(v))==1)))
        stop(sprintf("'%s' should be a character of length 1.", v))

    for(v in c("force", "do.logtransform", "grouprep", "spatial"))
      if (!(is.logical(get(v)) && (length(get(v))==1)))
        stop(sprintf("'%s' should be a logical of length 1.", v))

    if(!(intgroup %in% colnames(pData(expressionset))))
      stop("'intgroup' should match one of the column names of 'phenoData(expressionset)'.")


    ## Get going:
    olddir = getwd()
    on.exit(setwd(olddir))
    dircreation(outdir, force)
                
    arg = as.list(match.call(expand.dots = TRUE))
    name = arg$expressionset
    name = deparse(substitute(name))
    obj = list()
    
    datap = aqm.prepdata(expressionset, do.logtransform, sN)

    obj$maplot = try(aqm.maplot(dataprep = datap))
    if(inherits(obj$maplot, "try-error"))
      warning("Could not draw MA plots \n")

    if(inherits(expressionset, 'BeadLevelList') || inherits(expressionset, 'AffyBatch') ||
       ("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata)) && spatial == TRUE) {            
      obj$spatial =  try(aqm.spatial(expressionset = expressionset, dataprep = datap, scale = "Rank"))
      if(inherits(obj$spatial,"try-error"))
        warning("Could not draw spatial distribution of intensities \n")
                
      if((inherits(expressionset,'NChannelSet') && ("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset))))
         || inherits(expressionset,'BeadLevelList')) {
        obj$spatialbg =  try(aqm.spatialbg(expressionset = expressionset, dataprep = datap, scale = "Rank"))
        if(inherits(obj$spatialbg,"try-error"))
          warning("Could not draw spatial distribution of background intensities \n") 
      }                
    }
                        
    obj$boxplot = try(aqm.boxplot(expressionset = expressionset, dataprep = datap, intgroup, grouprep))
    if(inherits(obj$boxplot,"try-error"))
      warning("Could not draw boxplots \n")

    obj$density = try(aqm.density(expressionset = expressionset, dataprep = datap, intgroup, grouprep, outliers = unlist(obj$boxplot$outliers)))
    if(inherits(obj$density,"try-error"))
      warning("Could not draw density plots \n")

    obj$heatmap = try(aqm.heatmap(expressionset = expressionset, dataprep = datap, intgroup))
    if(inherits(obj$heatmap,"try-error"))
      warning("Could not draw heatmap \n") 

    if(all(intgroup %in% names(phenoData(expressionset)@data))) {
      obj$pca = try(aqm.pca(expressionset = expressionset, dataprep = datap, intgroup))
      if(inherits(obj$pca,"try-error"))
        warning("Could not draw PCA \n")
    }

    obj$meansd = try(aqm.meansd(dataprep = datap))
    if(inherits(obj$meansd,"try-error"))
      warning("Could not draw Mean vs Standard Deviation \n")
          

    if(inherits(expressionset,'BeadLevelList'))
      warning("Could not plot the probes mapping densities on a BeadLevelList object.")
    if(!inherits(expressionset,'BeadLevelList')) {
      if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
        {
          obj$probesmap = try(aqm.probesmap(expressionset = expressionset, dataprep = datap))
          if(inherits(obj$probesmap,"try-error"))
            warning("Cannot draw probes mapping plot \n")
        }
    }
    
    
    if(inherits(expressionset, "AffyBatch"))
      {
        obj$rnadeg = try(aqm.rnadeg(expressionset, sN))
        if(inherits(obj$rnadeg,"try-error"))
          warning("Cannot draw the RNA degradation plot \n")
        
        affyproc = aqm.prepaffy(expressionset, datap$sN)
        obj$rle = try(aqm.rle(affyproc))
        if(inherits(obj$rle,"try-error"))
          warning("Cannot draw the RLE plot \n") 
        
        obj$nuse = try(aqm.nuse(affyproc))
        if(inherits(obj$nuse,"try-error"))
          warning("Cannot draw the NUSE plot \n") 
        
        if(length(grep("exon", cdfName(expressionset), ignore.case=TRUE)) == 0)
          {
            obj$qcstats =  try(aqm.qcstats(expressionset))
            if(inherits(obj$qcstats,"try-error"))
              warning("Cannot draw the QCStats plot \n")
          }
        
        obj$pmmm = try(aqm.pmmm(expressionset))
        if(inherits(obj$pmmm,"try-error"))
          warning("Cannot draw the Perfect Match versus MisMatch plot \n") 
      }
    
    for(i in seq(along = obj))
      obj[[i]]$legend = gsub("The figure <!-- FIG -->",paste("<b>Figure",i, "</b>"),obj[[i]]$legend, ignore.case = TRUE) 
    
    aqm.writereport(name, expressionset, obj)
    return(invisible(obj))
  })



## RGList
setMethod("arrayQualityMetrics",signature(expressionset = "RGList"), function(expressionset, outdir, force, do.logtransform, intgroup, grouprep, spatial, sN)
          {
            expressionset = try(as(expressionset, "NChannelSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a RGList and cannot be converted automatically in a NChannelSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, intgroup, grouprep, spatial, sN)
          })####end set method RGList

## marrayRaw
setMethod("arrayQualityMetrics",signature(expressionset = "marrayRaw"), function(expressionset, outdir, force, do.logtransform, intgroup, grouprep, spatial, sN)
          {
            expressionset = try(as(expressionset, "NChannelSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a marrayRaw and cannot be converted automatically in a NChannelSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, intgroup, grouprep, spatial, sN)
          })####end set method marrayRaw

## MAlist
setMethod("arrayQualityMetrics",signature(expressionset = "MAList"), function(expressionset, outdir, force, do.logtransform, intgroup, grouprep, spatial, sN)
          {
            expressionset = try(new("ExpressionSet", exprs = expressionset$M, phenoData=new("AnnotatedDataFrame",expressionset$targets), featureData=new("AnnotatedDataFrame",expressionset$genes)))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a MAList and cannot be converted automatically in a ExpressionSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, intgroup, grouprep, spatial, sN)
          })####end set method MAList

## marrayNorm
setMethod("arrayQualityMetrics",signature(expressionset = "marrayNorm"), function(expressionset, outdir, force, do.logtransform, intgroup, grouprep, spatial, sN)
          {
            expressionset = try(as(expressionset, "ExpressionSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a marrayNorm and cannot be converted automatically in a ExpressionSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, intgroup, grouprep, spatial, sN)
          })####end set method marrayNorm

