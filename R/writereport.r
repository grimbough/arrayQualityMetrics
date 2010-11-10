## Creation of the outdir
dircreation = function(outdir = getwd(), force = FALSE)
  {
    if(file.exists(outdir)){
      if(!file.info(outdir)$isdir)
        stop(sprintf("'%s' must be a directory.", outdir))
      
      outdirContents = dir(outdir, all.files = TRUE)
      outdirContents = setdiff(outdirContents, c(".", ".."))
              
      if(!force && length(outdirContents)>0) 
        stop(sprintf("'%s' is not empty.", outdir))
      message(sprintf("The report will be written into directory '%s'. ", outdir))
    } else {
      dir.create(outdir, recursive=TRUE)
      message(sprintf("The directory '%s' has been created.", outdir))
    }
  }


## Produce the plots 
aqm.plot = function(x) {
  if (is(x@plot, "trellis") || is(x@plot, "list")) ## TODO remove 'list' once maplot is fixed
    print(x@plot) else
  if (is(x@plot, "function"))
    do.call(x@plot, args = list())
  else stop(sprintf("Invalid 'x@plot' of class '%s'.", paste(class(x@plot), collapse=", ")))
}

## Create the title
aqm.make.title = function(reporttitle, outdir, params)
  {
    if(!(is.character(reporttitle)&&(length(reporttitle)==1)))
      stop("'reporttitle' must be a character of length 1.")

    ## ------- copy and link the CSS and JavaScript files
    filenames = c("arrayQualityMetrics.css", "arrayQualityMetrics.js")
    filelocs  = system.file("javascript", filenames, package = "arrayQualityMetrics")
    filelocs[ filelocs!="" ]
    if(length(filelocs)<length(filenames))
        stop(sprintf("Could not find all of: '%s'.", paste(filenames, collapse=", ")))
    
    copySubstitute(src = filelocs, dest = outdir, recursive = TRUE, symbolValues = params)
                
    p = aqm.hwriteOpenPage(file.path(outdir, 'QMreport.html'),
      link.javascript = filenames[2],
      link.css        = filenames[1],
      body.attributes = c("onload" = "reportinit()"))

    hwrite("<hr>", p)
    hwrite(reporttitle, p, heading=1)
    hwrite("<hr>", p)
    return(p)
  }


## Create a new section
aqm.make.section = function(p, s, qm)
  {
    hwrite("<hr>", p)
    sec = paste("<a name= 'S",s,"'>Section ", s, ": ", qm@section,"</a>", sep = "")
    hwrite(sec, p, heading=2)
  }

## Create the index
aqm.make.index = function(modules, p)
{
  s = 1
  hwrite("<hr>", p)
  hwrite("Index", p, heading=2, style='font-family:helvetica,arial,sans-serif')
  hwrite("Note: bitmap figures shown in this HTML report are also linked to PDF files, to provide better figure quality, or to provide multi-page files in cases where in  this HTML report only the first page is presented.", p, style='font-family:helvetica;font-size:11pt;color:#808080;font-style:italic;')

  hwrite("<UL>", p)
  lasttype = "FAKE"
  for(i in seq_len(length(modules)))
    {
      if(modules[[i]]@section != lasttype)
        {
          if(s != 1) ## end the previous section
            hwrite("</UL>", p)

          hwrite(paste("<br><LI>Section ", s,": ", modules[[i]]@section, "<UL>",sep=""), p,
                 link = paste("#S",s,sep=""),
                 style = 'font-weight:bold;font-family:helvetica;font-size:12pt')
          s = s+1
        }
      hwrite(paste("<LI>", modules[[i]]@title,sep=""), p,
             style = 'font-weight:normal;font-family:helvetica;font-size:11pt')
      lasttype = modules[[i]]@section
    }
  hwrite("</UL></UL>", p)
}

##--------------------------------------------------
## 'tooltips' table for the mouseover events
##--------------------------------------------------
annotationTable = function(aTab, id, width=300) {
  name = paste(id, "cell", sep=":")
  mat  = cbind(paste("<td>", c("rowname", colnames(aTab)), "</td>", sep=""), "<td name='", name,"'></td>")
  rows = paste("<tr>", apply(mat, 1, paste, collapse=""), "</tr>", sep="")
  tab  = paste(rows, collapse="\n")
  tab  = paste("<table id='", id, "' width=", width, ">", tab, "</table>", sep="") 
}

## Create a part of report with figures and legend
aqm.report.qm = function(p, qm, f, name, arrayTable, outdir)
  {
    stopifnot(is(qm, "aqmReportModule"))

    h = qm@shape$h
    w = qm@shape$w
    dpi = 72

    imageformat =  if(length(qm@svg)>0) "svg" else "png"
    svgwarn = FALSE

    figID = paste("Fig", name, sep=":")
    tabID = paste("Tab", name, sep=":")
    
    ## The two cases, png and svg, need to be treated differently 
    switch(imageformat,
           svg = {
             nameimg = paste(name, ".svg", sep = "")
             svgtemp = tempfile()
             svg(file = svgtemp, h = h, w = w)
             aqm.plot(qm)
             dev.off()

              -------- Here we need to continue cleaning up -----
             annRes = annotateSvgPlot(infile=svgtemp, outfile=nameimg, outdir=outdir,
               annotationInfo = c(qm@svg, tabID = tabID)) 
             if(!annRes$annotateOK)
               svgwarn = "Note: the figure is static - enhancement with interactive effects failed. This is likely due to a version incompatibility of the arrayQualityMetrics package and libcairo. Please contact the Bioconductor mailing list to report this problem." 
             sizes = paste(annRes$size)
             img = hwrite(c(aqm.hwriteImage(nameimg, width=sizes[1], height=sizes[2], id=figID),
                            annotationTable(arrayTable, id=tabID)))
           },
           png = {
             nameimg = paste(name, ".png", sep = "")
             png(file = file.path(outdir, nameimg), h = h*dpi, w = w*dpi)
             aqm.plot(qm)
             dev.off()
             img = aqm.hwriteImage(nameimg, id=figID)
           },
           stop(sprintf("Invalid value of 'imageformat': %s", imageformat))
           )
    
    namepdf = paste(name, ".pdf", sep = "")
    pdf(file = file.path(outdir, namepdf), h = h, w = w)
    aqm.plot(qm)
    dev.off()
    link = list(namepdf, NA)
      
    linkpdf = hwrite("PDF file.",
      style='font-weight:bold;text-align:center;font-family:helvetica',
      border = 0, link = namepdf)
    
    hwrite(c(img, paste("Figure ", f, ": ", qm@title,". ", linkpdf, sep="")),
           p,
           dim=c(2,1),
           style='font-weight:bold;text-align:center;font-family:helvetica',
           border=0)


    hwrite(paste("<br>", gsub("The figure <!-- FIG -->", paste("<b>Figure", f, "</b>"), qm@legend, ignore.case = TRUE)), 
           p,
           style='text-align:justify;font-family:Lucida Grande;font-size:10pt')

    if(!identical(svgwarn, FALSE))
       hwrite(svgwarn,
           p,
           style='text-align:justify;font-family:Lucida Grande;font-size:10pt;color:blue')
    
  }

## End the report
aqm.make.ending = function(p)
  {
    z = sessionInfo("arrayQualityMetrics")
    version = z$otherPkgs[[1]]$Version
    rversion = sessionInfo()$R.version$version.string
    session = paste("This report has been created with arrayQualityMetrics", version, "under", rversion)
    hwrite("<hr>", p)
    hwrite(session, p, style ='font-family:Lucida Grande;font-size:8pt')
    hwrite("<hr>", p)
    closePage(p)
  }

## Score table formatting
scores = function(obj)
  {
    titles = c(
      "maplot"  = "MA plots",
      "spatial" = "Spatial",
      "boxplot" = "Distribution",
      "heatmap" = "Heatmap",
      "rle"     = "RLE",
      "nuse"    = "NUSE")

    mt = match(names(obj), names(titles))

    df = data.frame(
      "Array #"    = seq_len(length(sN)),
      "Array Name" = sN,
      stringsAsFactors = FALSE)

    ## TODO - this should go into vignette, explain how to put this into the pData.
    if(!is.null(protocolData(expressionset)$ScanDate))
      df$"Scan Dates" = protocolData(expressionset)$ScanDate
    
    for(i in which(!is.na(mt)))      
      df[[titles[mt[i]]]] =
        ifelse(seq(along=sN) %in% obj[[i]]@outliers, "*", "")

    return(df)
  }

aqm.make.table = function(arrayTable, p) {

   ## TO DO: add outlier detection
  
  arrayTable = cbind(
    highlight = sprintf("<input type=\"checkbox\" name=\"ArrayCheckBoxes\" value=\"aqm_%d\" onchange=\"checkboxEvent(%d)\"/>",
                        0:(nrow(arrayTable)-1),  0:(nrow(arrayTable)-1)),
    arrayTable, stringsAsFactors = FALSE)
  
  hwrite(arrayTable, p, border=0,
         row.bgcolor = rep(list("#ffffff", c("#d0d0ff", "#e0e0f0")), ceiling(nrow(arrayTable)/2)),
         cellpadding = 2, cellspacing = 5,
         row.style = list('font-weight:bold'))

  hwrite("<div id=\"arraytooltip\"></div>", p)

  
}


##--------------------------------------------------
##   write the report
##--------------------------------------------------

aqm.writereport = function(modules, arrayTable, reporttitle, outdir)
  {
    ## To avoid dealing with this pathologic special case downstream in the HTML
    if(nrow(arrayTable)==0)
      stop("'arrayTable' must not be empty.")

    ## For all report modules, extra the 'svg' slot, then subset only those that are defined.
    svgdata = lapply(modules, slot, "svg")
    svgdata = svgdata[ sapply(svgdata, slot, "defined") ]

    ## Extract strokeopacity and strokewidth from the list 'svgdata' and format for Javascript
    ##   (second part could also be done by RJSONIO)
    formatStrokeParameters = function(w)
      paste("[", apply(sapply(svgdata, slot, w), 2, paste, collapse=", "), "]", collapse=", ")

    ## Could also use RJSONIO here
    df = cbind( " " = rownames(arrayTable), arrayTable, stringsAsFactors=FALSE)
    pDataJS = sapply(df, function(x) paste("'", x, "'", sep=""))
    pDataJS = paste("[", apply(pDataJS, 1, paste, collapse=", "), "]", sep="", collapse=", ")

    ## Open and set up the HTML page
    ## Inject report-specific variables into the JavaScript
    p = aqm.make.title(
      reporttitle = reporttitle,
      outdir = outdir,
      params = c(          ## TODO: this should be the outliers
        HIGHLIGHTINITIAL = paste(rep(c("false", "true"), ceiling(nrow(arrayTable)/2))[seq_len(nrow(arrayTable))], collapse=", "),
        ARRAYMETADATA    = pDataJS,
        SVGOBJECTIDS     = paste("'Fig:", names(svgdata), "'", sep="", collapse=", "),
        TABLEIDS         = paste("'Tab:", names(svgdata), "'", sep="", collapse=", "),
        IDFUNS           = paste(sapply(svgdata, slot, "idFun"), collapse=", "),
        STROKEOPACITY    = formatStrokeParameters("strokeopacity"),
        STROKEWIDTH      = formatStrokeParameters("strokewidth")))

    aqm.make.table(arrayTable, p)
    
    aqm.make.index(modules, p)
    lasttype = "Something Else"
    sec = 1

    for(i in seq(along = modules))
      {
        if(modules[[i]]@section != lasttype)
          {
            aqm.make.section(p, s = sec, qm = modules[[i]])
            sec = sec+1
          }
        aqm.report.qm(p, qm=modules[[i]], f=i, name=names(modules)[i], arrayTable=arrayTable, outdir=outdir)
        lasttype = modules[[i]]@section
      }
     
    aqm.make.ending(p)
                    
    ## TODO is there a more elegant way to do this?
    invisible(list(modules=modules, arrayTable=arrayTable, reporttitle=reporttitle, outdir=outdir))
  }
