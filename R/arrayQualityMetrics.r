arrayQualityMetrics = function(
  expressionset,
  outdir = getwd(),
  force = FALSE,
  do.logtransform = FALSE,
  intgroup = NULL,
  grouprep,
  spatial = TRUE,
  reporttitle = paste("Quality metrics report for", deparse(substitute(expressionset))),
  usesvg) {
    
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

  if(missing(usesvg)){
    ## Note: assignment within the if-condition
    if(! (usesvg <- capabilities()["cairo"]) )
      warning("capabilities()[\"cairo\"] is FALSE - all graphics will be static. Please install the cairo library for your R to obtain interactive svg graphics.") 
  } else {
    if( is.logical(usesvg) && (length(usesvg)==1) && !is.na(usesvg) )
      stop("'usesvg' must be TRUE or FALSE")
    if(!capabilities()["cairo"])
      stop("capabilities()[\"cairo\"] is FALSE - cannot produce interactive svg graphics. Please install the cairo library for your R.") 
  }

  ## output directory
  dircreation(outdir, force)
                
  obj = list()
  
  ## create a comprehensive data object 'x', with the original data,
  ##  as well as some generally useful derived statistics of the data
  x = prepdata(expressionset, intgroup, do.logtransform) 
  
  ## TODO: move this into function 'aqm.spatial'
  if(inherits(expressionset, 'BeadLevelList') || inherits(expressionset, 'AffyBatch') ||
     (("X" %in% rownames(featureData(expressionset)@varMetadata)) &&
      ("Y" %in% rownames(featureData(expressionset)@varMetadata))) && spatial) {            
    obj$spatial =  try(aqm.spatial(expressionset = expressionset, dataprep = x, scale = "Rank"))
    if(inherits(obj$spatial,"try-error"))
      warning("Could not draw spatial distribution of intensities \n")
    
    if((inherits(expressionset,'NChannelSet') &&
        ("Rb" %in% colnames(dims(expressionset)) &&
         ("Gb" %in% colnames(dims(expressionset)))))
       || inherits(expressionset,'BeadLevelList')) {
      obj$spatialbg =  try(aqm.spatialbg(expressionset = expressionset, dataprep = x, scale = "Rank"))
      if(inherits(obj$spatialbg,"try-error"))
        warning("Could not draw spatial distribution of background intensities \n") 
    }                
  }
  
  obj$maplot    = aqm.maplot (x)
  obj$boxplot   = aqm.boxplot(x)
  obj$density   = aqm.density(x, outliers = obj$boxplot$outliers)
  obj$heatmap   = aqm.heatmap(x)
  obj$pca       = aqm.pca    (x, outliers = obj$heatmap$outliers)
  obj$meansd    = aqm.meansd (x)
  obj$probesmap = aqm.probesmap(x)

  ##--------Affymetrix specific things------------
  if(inherits(expressionset, "AffyBatch")) {
    affyproc   = prepaffy(expressionset)
    obj$rnadeg = aqm.rnadeg(expressionset)
    obj$rle    = aqm.rle(x,  affyproc = affyproc)
    obj$nuse   = aqm.nuse(x, affyproc = affyproc)
    if(length(grep("exon", cdfName(expressionset), ignore.case=TRUE)) == 0)
      obj$qcstats =  aqm.qcstats(expressionset)
    
    obj$pmmm = aqm.pmmm(expressionset)
  }
  
  remove = names(obj)[ sapply(obj, inherits, what = "try-error") ]
  for(r in remove)
    obj[[r]] = NULL
  
  for(i in seq(along = obj))
    obj[[i]]$legend = gsub("The figure <!-- FIG -->", paste("<b>Figure", i, "</b>"), obj[[i]]$legend, ignore.case = TRUE)
  
  aqm.writereport(name = reporttitle, obj = obj, outdir = outdir)
  return(invisible(obj))
}

