## --------------------------------------------------------------------------------------------------------
## Postprocess an SVG file:
##  1. add an 'id' attribute to the root svg element
##  2. the ids of the symbols that are defined in the <defs> elements are -by libcairo- of the form
##     glyph0-0, glyph1-0 etc. Since we are going to inline the content from all svg files into a single html
##     document, these id would clash between different svg plots, and need to by made unique. I wonder
##     whether a more elegant way exists for this.
##  3. add mouse events to the elements of interest (found by the function annotationInfo@getPlotObjNodes)
## ---------------------------------------------------------------------------------------------------------
annotateSvgPlot = function(infile, outfile, outdir, annotationInfo, name) 
  {
     
    ## Check argument 
    stopifnot(is(annotationInfo, "svgParameters"))
  
    doc = xmlParse(infile)
    vb  = getViewBox(doc)

    svg = xmlRoot(doc)
    
    ## 1. add id
    xmlAttrs(svg) = c(id = paste("Fig", name, sep=":"))  

    ## monitor our success in finding what we expect
    isok = c(syms = FALSE, uses = FALSE, plotobjs = FALSE)
    
    ## 2. find the children of the <defs> element that are <symbol> 
    syms = getNodeSet(doc, "//x:defs/x:g/x:symbol", "x") 
    if(length(syms)>0)
      {
        oldids = sapply(syms, function(x) xmlAttrs(x)["id"])
        if(all(grepl("^glyph", oldids)))
          {
            newids = sub("^glyph", paste(name, "-", sep=""), oldids)
            for(i in seq(along=syms))
              xmlAttrs(syms[[i]]) = newids[i] 
            isok["syms"] = TRUE
          }
      }
    
    ## similarly, find the <use> elements 
    uses = getNodeSet(doc, "//x:use", "x") 
    if(length(uses)>0)
      {
        oldhrefs = sapply(uses, function(x) xmlAttrs(x)["href"])
        if(all(grepl("^#glyph", oldhrefs)))
          {
            newhrefs = sub("^#glyph", paste("#", name, "-", sep=""), oldhrefs)
            names(newhrefs) = rep("xlink:href", length(newhrefs))
            for(i in seq(along=uses))
              xmlAttrs(uses[[i]]) = newhrefs[i]
            isok["uses"] = TRUE
          }
      }

    browser()
    
    ## 3. this part is brittle - 'getPlotObjNodes' will be 'getMatplotSeries' or 'getPlotPoints' from
    ## 'SVGAnnotation', which rely on conventions used by libcairo to produce the SVG
    ## from the R plot, on simple pattern matching and on hope that the found patterns
    ## align with the intended plot objects (i.e. not on any explicit identification).
    plotobjs = annotationInfo@getPlotObjNodes(doc)
    
    if (length(plotobjs) == annotationInfo@numPlotObjects)
      {
        for(i in seq(along=plotobjs))
          {
            poid = paste("p", i, sep=":")
            roid = annotationInfo@getReportObjIdFromPlotObjId(poid)
            stopifnot(length(roid)==1)
            
            callbacks = sprintf("top.plotObjRespond('%s', '%s', '%s', '%s')", c("click", "show", "hide"), poid, roid, name)
            
            xmlAttrs(plotobjs[[i]]) = c(
                      "id"          = poid,
                      "onclick"     = callbacks[1],
                      "onmouseover" = callbacks[2],
                      "onmouseout"  = callbacks[3])
            
            convertCSSStylesToSVG(plotobjs[[i]])
          } ## for
        isok["plotobjs"] = TRUE
      }
    
    saveXML(doc, file.path(outdir, outfile))
    return(list(size = diff(vb), annotateOK = all(isok)))
  }

##--------------------------------------------------------------------------------------
## HTML table to show 'tooltips' for mouseover events
## The function creates a table with 2 columns and as many rows as 'x' has columns.
## The first column will contain the rownames of 'x', the second column will be empty
##---------------------------------------------------------------------------------------
annotationTable = function(x, name) {
  bgcol = rep(c("#d0d0ff", "#e0e0f0"), ceiling(ncol(x)/2))[seq_len(ncol(x))]
  tab  = paste("<tr bgcolor='", bgcol, "'><td>", colnames(x), "</td><td style='font-weight:bold'></td></tr>\n", sep="", collapse="\n")
  tab  = paste("<table id='", paste("Tab", name, sep=":"), "'>", tab, "</table>", sep="")
}

