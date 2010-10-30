## postprocess svg file to provide explanation on objects
##
## annotation: named list, length matches the number of elements to be annotated
##   the names are the unique ids
##   each list element is again a list of
##     - title (character string)
##     - linked ids (character vector)
##

annotateSvgMatplot = function(infile, outfile, outdir, annotationInfo) 
  {
     
    doc = xmlParse(infile)
    svg = xmlRoot(doc)
    vb  = getViewBox(doc)
    
    ## Add the lineWidth toggling
    annotateOK = aqm.highlight(doc, annotationInfo = annotationInfo)
    
    ## Did it work?
    if(annotateOK){
      
      ## Enlarge SVG view box to make spake for the "status bar"
      addy = 20
      enlargeSVGViewBox(doc, 0, addy) ## nothing in x-direction, 100 units in y-direction
      
      ## Add the "status bar"
      vb = getViewBox(doc)
      g = xmlRoot(doc)[["g"]]
      newXMLNode("text", " ",
                 attrs = c("x" = vb[1,1],
                   "y" = vb[2,2]-addy,
                   "font-size" = "13",
                   "font-family" = "Arial,Helvetica",
                   "id" = "annotationtext"),
                 parent = g)
      
      ## addCSS(doc, insert = TRUE)

      oldwd = setwd(outdir)  ## the following function needs this
      addECMAScripts(doc, scripts = "arrayQualityMetrics.js", insertJS = FALSE)
      setwd(oldwd)
      
      ## Add an onload Java script call to the <svg> tag
      addAttributes(svg, "onload"="init(evt);")
    }
    
    saveXML(doc, file.path(outdir, outfile))

    return(list(size=diff(vb), annotateOK=annotateOK))
  }

## The following is adapted from the functions
## highlightMatplot and highlightMatplotSeries in the package
## SVGAnnotation
aqm.highlight = function(doc, annotationInfo)
{
  
  series = annotationInfo$getfun(doc)
  anno = annotationInfo$annotation
  
  ## TODO - this should never happen
  if(length(anno) != length(series)) {
    ## browser()
    return(FALSE)
  }
  
  for(i in seq(along=series)){
    ops = sprintf("toggleSeries(%s, %s, %s)",
      paste("[", paste("'", anno[[i]]$linkedids, "'", sep="", collapse=","), "]", sep=""),
      paste("'", anno[[i]]$title, "'", sep=""),
      c("true", "false"))
    names(ops) = c("onmouseover", "onmouseout")
    
    node = series[[i]]
    xmlAttrs(node) = c(id = names(anno)[i], ops)
    convertCSSStylesToSVG(node)
    series[[i]] = node
  }
  return(TRUE)
}

## based upon SVGAnnotation::getMatplotSeries
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

  
