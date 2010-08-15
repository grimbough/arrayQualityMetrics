## postprocess svg file to provide explanation on objects

annotateSvgMatplot = function(infile,
  outfile,
  name,
  js = system.file("javaScript", "imatplot.js", package = "arrayQualityMetrics"),
  js = "imatplot.js") {
  
  doc = xmlParse(infile)
  svg = xmlRoot(doc)

  addECMAScripts(doc, js, insertJS = TRUE)
  addCSS(doc, insert = TRUE)

  ## Add the lineWidth toggling
  aqm.highlightMatplot(doc, name = name)

  ## Enlarge SVG view box to make spake for the "status bar"
  addy = 20
  enlargeSVGViewBox(doc, 0, addy) ## nothing in x-direction, 100 units in y-direction

  ## Add the "status bar"
  vb = getViewBox(doc)
  g = xmlRoot(doc)[["g"]]
  newXMLNode("text", "Horizontal",
              attrs = c("x" = vb[1,1]+50,
                        "y" = vb[2,2]-addy,
                        "font-size" = "14",
                        "font-family" = "Arial,Helvetica",
                        "id" = "annotationtext"),
              parent = g)
  
  ## Add an onload Java script call to the <svg> tag
  addAttributes(svg, "onload"="init(evt);")
  
  saveXML(doc, outfile)
}

## The following two functions are adapted from the functions
## highlightMatplot and highlightMatplotSeries in the package
## SVGAnnotation
aqm.highlightMatplot= function(doc,
    series = getMatplotSeries(doc),
    id)
{
    if(length(id) != length(series))
      stop("'id' must be equal to 'length(series)' (the number of lines in the plot).")
    mapply(aqm.highlightMatplotSeries, series, id = id)
}

aqm.highlightMatplotSeries = function (node, id) 
{
    ops = sprintf("toggleSeries(\"%s\", %s)", id, c("true", "false"))
    names(ops) = c("onmouseover", "onmouseout")
    xmlAttrs(node) = c(id = id, ops)
    convertCSSStylesToSVG(node)
}
