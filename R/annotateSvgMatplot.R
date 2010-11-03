## --------------------------------------------------
## Postprocess an SVG file to add mouse events
## --------------------------------------------------
annotateSvgPlot = function(infile, outfile, outdir, annotationInfo) 
  {
     
    doc = xmlParse(infile)
    svg = xmlRoot(doc)
    vb  = getViewBox(doc)
    
    annotateOK = annotatePlotObjects(doc, annotationInfo)
    
    saveXML(doc, file.path(outdir, outfile))
    return(list(size=diff(vb), annotateOK = annotateOK))
  }

## --------------------------------------------------------------------------
## Find the individual objects (lines, points) in the plot and add the events
## --------------------------------------------------------------------------
annotatePlotObjects = function(doc, annotationInfo)
{
  ## extract and check arguments 
  getfun = annotationInfo$getfun
  numObj = annotationInfo$numObjects
  tabID  = annotationInfo$tabID
  stopifnot( !is.null(getfun), is.function(getfun),
             !is.null(numObj), is.numeric(numObj),  length(numObj)==1, !is.na(numObj),
             !is.null(tabID),  is.character(tabID), length(tabID)==1,  !is.na(tabID) )
  
  ## This part is brittle - 'getfun' will be 'getMatplotSeries' or 'getPlotPoints' from
  ## SVGAnnotation, which rely on conventions used by libcairo to produce the SVG
  ## from the R plot, on simple pattern matching and on hope that the found patterns
  ## align with the intended plot objects (i.e. not on any explicit identification).
  series = annotationInfo$getfun(doc)
  
  ## Catch some of the brittleness
  if( (length(series) %% numObj) != 0) {
    return(FALSE)
  }

  if (length(series)>0)
    for(pel in 0:(length(series)-1))
      {
        ## Use integer 'pel' to refer to the plot elements (line or symbol; in this package: corresponding to an array).
        ops = c("onclick" = sprintf("top.clickPlotElement(%d)", pel %% numObj),
            "onmouseover" = sprintf("top.showTip('%s', %d)", tabID, pel %% numObj),
             "onmouseout" = sprintf("top.hideTip('%s')", tabID ))
        
        node = series[[pel+1]]     ## some pain dealing with 0 versus 1 based arrays...
        xmlAttrs(node) = c(id = sprintf("aqm%d", pel), ops)
        convertCSSStylesToSVG(node)
        series[[pel+1]] = node
      }
  return(TRUE)
}

##-----------------------------------------------------------------
## Check the 'usesvg' parameter, and the available infrastructure,
## and if appropriate, emit a warning.
##-----------------------------------------------------------------
checkUsesvg = function(usesvg) {

  if(missing(usesvg)){
    ## Note: assignment within the if-condition
    if(! (usesvg <- capabilities()["cairo"]) )
      warning("capabilities()[\"cairo\"] is FALSE - all graphics will be static. Please install the cairo library for your R to obtain interactive SVG graphics.") 
  } else {
    if( is.logical(usesvg) && (length(usesvg)==1) && !is.na(usesvg) )
      stop("'usesvg' must be TRUE or FALSE")
    if(usesvg && (!capabilities()["cairo"]))
      stop("capabilities()[\"cairo\"] is FALSE - cannot produce interactive SVG graphics. Please install the cairo library for your R.") 
  }
  
  return(usesvg)
}

  
