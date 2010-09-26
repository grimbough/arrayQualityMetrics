setClassUnion("aqmInputObj", c("ExpressionSet", "AffyBatch", "NChannelSet", "BeadLevelList"))

## set the method for arrayQualityMetrics
setGeneric("arrayQualityMetrics",
  def = function(expressionset,
           outdir = getwd(),
           force = FALSE,
           do.logtransform = FALSE,
           intgroup = NULL,
           grouprep,
	   spatial = TRUE,
           reporttitle = paste("Quality metrics report for", deparse(substitute(expressionset)))
  )
  standardGeneric("arrayQualityMetrics"),
  signature = "expressionset")

## aqmInputObj
setMethod("arrayQualityMetrics",
          signature(expressionset = "aqmInputObj"),
          function(expressionset,
                   outdir = getwd(),
                   force = FALSE,
                   do.logtransform = FALSE,
                   intgroup = NULL,
                   grouprep,
                   spatial = TRUE,
                   reporttitle = paste("Quality metrics report for", deparse(substitute(expressionset)))) {
    
    ## Argument checking: 
    if(!missing(grouprep))
      .Deprecated(msg = paste("The argument 'grouprep' of the function 'arrayQualityMetrics'",
        "is deprecated and will be ignored. Use 'intgroup' instead."))
      
    for(v in c("outdir", "reporttitle"))
      if (!(is.character(get(v)) && (length(get(v))==1)))
        stop(sprintf("'%s' should be a character of length 1.", v))

    for(v in c("force", "do.logtransform", "spatial"))
      if (!(is.logical(get(v)) && (length(get(v))==1)))
        stop(sprintf("'%s' should be a logical of length 1.", v))

    if (!is.null(intgroup)) {
      if (!is.character(intgroup))
        stop("'intgroup' should be a 'character'.")
      if(!all(intgroup %in% colnames(pData(expressionset))))
        stop("all elements of 'intgroup' should match column names of 'phenoData(expressionset)'.")
    }

    ## Get going:
    olddir = getwd()
    on.exit(setwd(olddir))
    dircreation(outdir, force)
                
    obj = list()
    
    dataprep = aqm.prepdata(expressionset, do.logtransform)

    obj$maplot = aqm.maplot(dataprep = dataprep)

    if(inherits(expressionset, 'BeadLevelList') || inherits(expressionset, 'AffyBatch') ||
       (("X" %in% rownames(featureData(expressionset)@varMetadata)) &&
        ("Y" %in% rownames(featureData(expressionset)@varMetadata))) && spatial) {            
      obj$spatial =  try(aqm.spatial(expressionset = expressionset, dataprep = dataprep, scale = "Rank"))
      if(inherits(obj$spatial,"try-error"))
        warning("Could not draw spatial distribution of intensities \n")
                
      if((inherits(expressionset,'NChannelSet') &&
          ("Rb" %in% colnames(dims(expressionset)) &&
          ("Gb" %in% colnames(dims(expressionset)))))
         || inherits(expressionset,'BeadLevelList')) {
        obj$spatialbg =  try(aqm.spatialbg(expressionset = expressionset, dataprep = dataprep, scale = "Rank"))
        if(inherits(obj$spatialbg,"try-error"))
          warning("Could not draw spatial distribution of background intensities \n") 
      }                
    }
                        
    obj$boxplot = aqm.boxplot(expressionset, dataprep=dataprep, intgroup = intgroup)
    obj$density = aqm.density(expressionset, dataprep=dataprep, intgroup = intgroup, outliers = obj$boxplot$outliers)

    obj$heatmap = aqm.heatmap(expressionset, dataprep=dataprep, intgroup = intgroup)
    obj$pca     = aqm.pca(    expressionset, dataprep=dataprep, intgroup = intgroup, outliers = obj$heatmap$outliers)
    obj$meansd  = aqm.meansd( dataprep=dataprep)

    if(inherits(expressionset,'BeadLevelList'))
      warning("Could not plot the probes mapping densities on a BeadLevelList object.")
    else if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
          {
            obj$probesmap = try(aqm.probesmap(expressionset = expressionset, dataprep = dataprep))
            if(inherits(obj$probesmap,"try-error"))
              warning("Could not draw probes mapping plot \n")
          }
    
    
    if(inherits(expressionset, "AffyBatch"))
      {
        obj$rnadeg = try(aqm.rnadeg(expressionset))
        if(inherits(obj$rnadeg,"try-error"))
          warning("Could not draw the RNA degradation plot \n")
        
        affyproc = aqm.prepaffy(expressionset)
        obj$rle = try(aqm.rle(expressionset, dataprep, affyproc, intgroup=intgroup))
        if(inherits(obj$rle,"try-error"))
          warning("Could not draw the RLE plot \n") 
        
        obj$nuse = try(aqm.nuse(expressionset, dataprep, affyproc, intgroup=intgroup))
        if(inherits(obj$nuse,"try-error"))
          warning("Could not draw the NUSE plot \n") 
        
        if(length(grep("exon", cdfName(expressionset), ignore.case=TRUE)) == 0)
          {
            obj$qcstats =  try(aqm.qcstats(expressionset))
            if(inherits(obj$qcstats,"try-error"))
              warning("Could not draw the QCStats plot \n")
          }
        
        obj$pmmm = try(aqm.pmmm(expressionset))
        if(inherits(obj$pmmm,"try-error"))
          warning("Could not draw the Perfect Match versus Mismatch plot \n") 
      }
    
    remove = names(obj)[ sapply(obj, inherits, what = "try-error") ]
    for(r in remove) obj[[r]] = NULL
    
    for(i in seq(along = obj))
      obj[[i]]$legend = gsub("The figure <!-- FIG -->", paste("<b>Figure", i, "</b>"), obj[[i]]$legend, ignore.case = TRUE)

    aqm.writereport(reporttitle, expressionset, obj)
    return(invisible(obj))
  })



for(othertype in c("RGList", "MAList", "marrayRaw", "marrayNorm"))
  setMethod("arrayQualityMetrics", signature(expressionset = othertype),
            function(expressionset, outdir = getwd(), force = FALSE, do.logtransform = FALSE, intgroup = NULL, spatial = TRUE) {
              expressionset = try(as(expressionset, "NChannelSet"))
              if(inherits(expressionset,'try-error'))
                stop(sprintf("Argument 'expressionset' is of class '%s', and its automatic conversion into 'NChannelSet' failed. Please try to convert it manually.\n", othertype))
              arrayQualityMetrics(expressionset, outdir, force, do.logtransform, intgroup, spatial)
            })
