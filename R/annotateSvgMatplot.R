## --------------------------------------------------
## Postprocess an SVG file to add mouse events
## --------------------------------------------------
annotateSvgPlot = function(infile, outfile, outdir, annotationInfo, name) 
  {
     
    doc = xmlParse(infile)
    svg = xmlRoot(doc)
    vb  = getViewBox(doc)
    
    ## Extract and check arguments 
    stopifnot(is(annotationInfo, "svgParameters"))
  
    ## This part is brittle - 'getPlotObjNodes' will be 'getMatplotSeries' or 'getPlotPoints' from
    ## 'SVGAnnotation', which rely on conventions used by libcairo to produce the SVG
    ## from the R plot, on simple pattern matching and on hope that the found patterns
    ## align with the intended plot objects (i.e. not on any explicit identification).
    series = annotationInfo@getPlotObjNodes(doc)
  
    ## Try to catch some of the brittleness
    if (length(series) != annotationInfo@numPlotObjects)
      {
        annotateOK = FALSE
      }
    else
      {
        annotateOK = TRUE
        
        for(i in seq(along=series))
          {
            poid = paste("p", i, sep=":")
            roid = annotationInfo@getReportObjIdFromPlotObjId(poid)
            stopifnot(length(roid)==1)
            
            callbacks = sprintf("top.plotObjRespond('%s', '%s', '%s', '%s')", c("click", "show", "hide"), poid, roid, name)
            
            xmlAttrs(series[[i]]) = c(
                      "id"          = poid,
                      "onclick"     = callbacks[1],
                      "onmouseover" = callbacks[2],
                      "onmouseout"  = callbacks[3])
            
            convertCSSStylesToSVG(series[[i]])
          } ## for
      } ## else
    
    saveXML(doc, file.path(outdir, outfile))
    return(list(size = diff(vb), annotateOK = annotateOK))
  }

##--------------------------------------------------------------------------------------
## HTML table to show 'tooltips' for mouseover events
## The function creates a table with 2 columns and as many rows as 'x' has columns.
## The first column will contain the rownames of 'x', the second column will be empty
##---------------------------------------------------------------------------------------
## annotationTable = function(x, name, width=300) {
annotationTable = function(x, name) {
  bgcol = rep(c("#d0d0ff", "#e0e0f0"), ceiling(ncol(x)/2))[seq_len(ncol(x))]
  tab  = paste("<tr bgcolor='", bgcol, "'><td>", colnames(x), "</td><td style='font-weight:bold'></td></tr>", sep="", collapse="\n")
  tab  = paste("<table id='", paste("Tab", name, sep=":"), "'>", tab, "</table>", sep="")
  ## tab  = paste("<table id='", paste("Tab", name, sep=":"), "' width=", width, ">", tab, "</table>", sep="")
}

