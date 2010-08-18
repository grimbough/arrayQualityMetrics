## postprocess svg file to provide explanation on objects
##
## annotation: named list, length matches the number of elements to be annotated
##   the names are the unique ids
##   each list element is again a list of
##     - title (character string)
##     - linked ids (character vector)
##


annotateSvgMatplot = function(infile, outfile, annotationInfo,
  js = system.file("javascript", "imatplot.js", package = "arrayQualityMetrics")) 
  {
    stopifnot(js!="")
    
    doc = xmlParse(infile)
    svg = xmlRoot(doc)
    
    ## addCSS(doc, insert = TRUE)
    addECMAScripts(doc, js, insertJS = TRUE)
    
    ## Add the lineWidth toggling
    aqm.highlight(doc, annotationInfo = annotationInfo)
    
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
    
    ## Add an onload Java script call to the <svg> tag
    addAttributes(svg, "onload"="init(evt);")
    
    saveXML(doc, outfile)

    size = diff(vb)
    return(size)
  }

## The following is adapted from the functions
## highlightMatplot and highlightMatplotSeries in the package
## SVGAnnotation
aqm.highlight = function(doc, annotationInfo)
{
  
  series = annotationInfo$getfun(doc)
  anno = annotationInfo$annotation
  
  if(length(anno) != length(series))
    stop("'length(annotationInfo$annotation)' must be equal to 'length(series)', the number of objects in the plot.")
  
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
}
