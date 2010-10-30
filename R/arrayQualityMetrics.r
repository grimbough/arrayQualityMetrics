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

  usesvg = checkUsesvg( usesvg )
  
  ## output directory
  dircreation(outdir, force)

  ## list of report modules
  m = list()  
  
  ## create a comprehensive data object 'x', with the original data,
  ##  as well as some generally useful derived statistics of the data
  x = prepdata(expressionset, intgroup=intgroup, do.logtransform=do.logtransform, usesvg=usesvg) 
  
  ##---------Spatial intensity distributions------
  if (spatial) {
    hasxy = all(c("X", "Y") %in% rownames(featureData(expressionset)@varMetadata)) 
    if( inherits(expressionset, 'BeadLevelList') ||
        inherits(expressionset, 'AffyBatch') ||
        hasxy ) {
      
      m$spatial = aqm.spatial(expressionset = expressionset, dataprep = x, scale = "Rank")
    
      if((inherits(expressionset, 'NChannelSet') &&
        all(c("Gb"," Rb") %in% colnames(dims(expressionset)))) ||
         inherits(expressionset,'BeadLevelList')) {
        m$spatialbg =  aqm.spatialbg(expressionset = expressionset, dataprep = x, scale = "Rank")
      }
    }
  }
  
  ##---------Generic modules------
  m$maplot    = aqm.maplot (x)
  m$boxplot   = aqm.boxplot(x)
  m$density   = aqm.density(x, outliers = m$boxplot@outliers)
  m$heatmap   = aqm.heatmap(x)
  m$pca       = aqm.pca    (x, outliers = m$heatmap@outliers)
  m$meansd    = aqm.meansd (x)
  m$probesmap = aqm.probesmap(x)

  ##--------Affymetrix specific modules------------
  if(inherits(expressionset, "AffyBatch")) {
    m$rnadeg  = aqm.rnadeg(expressionset, x)
    affyproc  = prepaffy(expressionset)
    m$rle     = aqm.rle(x,  affyproc = affyproc)
    m$nuse    = aqm.nuse(x, affyproc = affyproc)
    if(length(grep("exon", cdfName(expressionset), ignore.case=TRUE)) == 0)
      m$qcstats =  aqm.qcstats(expressionset)
    
    m$pmmm = aqm.pmmm(expressionset)
  }
  
  for(i in seq(along = m))
    m[[i]]@legend = gsub("The figure <!-- FIG -->", paste("<b>Figure", i, "</b>"), m[[i]]@legend, ignore.case = TRUE)

  aqm.writereport(modules = m, arrayTable = x$pData, reporttitle = reporttitle, outdir = outdir)
}

