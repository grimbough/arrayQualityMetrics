arrayQualityMetrics = function(
  expressionset,
  outdir = reporttitle,
  force = FALSE,
  do.logtransform = FALSE,
  intgroup = NULL,
  grouprep,
  spatial = TRUE,
  reporttitle = paste("Quality metrics report for", deparse(substitute(expressionset))),
  usesvg)
{
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
  
  ##---------Generic modules------
  m$heatmap   = aqm.heatmap(x)
  m$pca       = aqm.pca    (x)

  m$boxplot   = aqm.boxplot(x)
  m$density   = aqm.density(x)
  
  m$meansd    = aqm.meansd (x)
  m$probesmap = aqm.probesmap(x)

  ##--------Affymetrix specific modules------------
  if(inherits(expressionset, "AffyBatch")) {
    x         = prepaffy(expressionset, x)
    m$rle     = aqm.rle(x)
    m$nuse    = aqm.nuse(x)
    m$rnadeg  = aqm.rnadeg(expressionset, x)
    ## m$qcstats = aqm.qcstats(expressionset)  -- Not sure any one cares about this function, anyway. It can be resurrected if there is overwhelming demand.
    m$pmmm    = aqm.pmmm(x)
  }
  
  ##---------MA plots and spatial intensity distributions------
  m$maplot = aqm.maplot (x)
  if (spatial) 
    m = append(m, aqm.spatial(x))
  
  aqm.writereport(modules = m, arrayTable = x$pData, reporttitle = reporttitle, outdir = outdir)
}

