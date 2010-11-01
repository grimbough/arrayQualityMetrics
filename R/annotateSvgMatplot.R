## --------------------------------------------------
## Postprocess an SVG file to add mouse events
## --------------------------------------------------
annotateSvgPlot = function(infile, outfile, outdir, annotationInfo) 
  {
     
    doc = xmlParse(infile)
    svg = xmlRoot(doc)
    vb  = getViewBox(doc)
    
    annotateOK = aqm.highlight(doc, annotationInfo)
    
    ## Did it work?
    if(annotateOK)
      {
        oldwd = setwd(outdir)  ## the following function needs this
        addECMAScripts(doc, scripts = "arrayQualityMetrics.js", insertJS = FALSE)
        setwd(oldwd)
      }
    saveXML(doc, file.path(outdir, outfile))
    return(list(size=diff(vb), annotateOK = annotateOK))
  }

## The following is adapted from the functions
## highlightMatplot and highlightMatplotSeries in the package
## SVGAnnotation
aqm.highlight = function(doc, annotationInfo)
{
  getfun = annotationInfo$getfun
  numObj = annotationInfo$numObjects
  
  stopifnot( !is.null(getfun), !is.null(numObj), is.numeric(numObj), length(numObj)==1, !is.na(numObj) )
  series = annotationInfo$getfun(doc)
  
  ## TODO - this should never happen
  if( (length(series) %% numObj) != 0) {
    ## browser()
    return(FALSE)
  }

  if (length(series)>0)
    for(i in 0:(length(series)-1))
      {
        ## For now, we use directly the integers 'i' to refer to the objects.
        ## This is OK since SVGAnnotation does not do anything more sophisticated anyway.
        ## Once the SVG files produced from R contain real identifiers for the plot objects (lines, points)
        ##   we should use these.
        ops = c("onclick" = sprintf("top.clickPlotElement(%d)", i %% numObj),
            "onmouseover" = sprintf("top.showTip(%d)", i %% numObj),
            "onmouseout"  = "top.hideTip()" )
        
        node = series[[i+1]]     ## some pain dealing with 0 versus 1 based arrays...
        xmlAttrs(node) = c(id = sprintf("aqm%d", i), ops)
        convertCSSStylesToSVG(node)
        series[[i+1]] = node
      }
  return(TRUE)
}

## This is one instance of a function that can be used for 'getfun'
## in aqm.highlight. It is based upon SVGAnnotation::getMatplotSeries.
## It used by aqm.density and aqm.rnadeg, which show "matplot"-like lines.
## aqm.pca, which shows a scatterplot with points, uses SVGAnnotation:::getPlotPoints
aqm.getMatplotSeries = function(doc)
    SVGAnnotation::getMatplotSeries(doc,
         paths = XML::getNodeSet(doc, "//x:g[starts-with(@id, 'surface')]//x:path", "x"))


## check the 'usesvg' parameter, and the available infrastructure, and if appropriate, emit warning.
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

  
