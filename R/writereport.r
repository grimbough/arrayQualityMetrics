## Creation of the output directory
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


## Produce a plots 
makePlot = function(x) {
  if (is(x@plot, "trellis") || is(x@plot, "list")) ## TODO remove 'list' once maplot is fixed
    print(x@plot) else
  if (is(x@plot, "function"))
    do.call(x@plot, args = list())
  else stop(sprintf("Invalid 'x@plot' of class '%s'.", paste(class(x@plot), collapse=", ")))
}

## Create the title
makeTitle = function(reporttitle, outdir, params)
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
                
    p = aqm.hwriteOpenPage(file.path(outdir, 'index.html'),
      link.javascript = filenames[2],
      link.css        = filenames[1],
      body.attributes = c("onload" = "reportinit()"))

    hwrite("<hr>", p)
    hwrite(reporttitle, p, heading=1)
    hwrite("<hr>", p)
    return(p)
  }


## Create a new section
makeSection = function(p, sectionNumber, module)
  {
    hwrite("<hr>", p)
    sec = paste("<a name= 'S", sectionNumber,"'>Section ", sectionNumber, ": ", module@section,"</a>", sep = "")
    hwrite(sec, p, heading=2)
  }

## Create the index
makeIndex = function(p, modules)
{
  s = 1
  hwrite("Index", p, heading=2, style='font-family:helvetica,arial,sans-serif')

  hwrite("<UL>", p)
  lasttype = "FAKE"
  for(i in seq_len(length(modules)))
    {
      if(modules[[i]]@section != lasttype)
        {
          if(s != 1) ## end the previous section
            hwrite("</UL>", p)

          hwrite(paste("<br><LI><div class='tocsection'>Section ", s,": ", modules[[i]]@section, "</div><UL>",sep=""), p,
                 link = paste("#S",s,sep=""))
          s = s+1
        }
      hwrite(paste("<LI><div class='tocmodule'>", modules[[i]]@title, "</div>", sep=""), p)
      lasttype = modules[[i]]@section
    }
  hwrite("</UL></UL>", p)
}

##---------------------------------------------------------
## Create a module  of the report with figures and legend
##---------------------------------------------------------
reportModule = function(p, module, integerIndex, name, arrayTable, outdir)
{
    stopifnot(is(module, "aqmReportModule"))

    h = module@shape$h
    w = module@shape$w
    stopifnot(length(h)==1, is.numeric(h), !is.na(h),
              length(w)==1, is.numeric(w), !is.na(w))

    stopifnot(length(integerIndex)==1, is.integer(integerIndex), !is.na(integerIndex))
    
    imageformat =  if(is.na(module@svg@name)) "png" else "svg"
    svgwarn = FALSE
    dpi = 72

    ## The two cases, png and svg, need to be treated differently 
    switch(imageformat,
        svg = {
          nameimg = paste(name, ".svg", sep = "")
          svgtemp = tempfile()
          svg(file = svgtemp, h = h, w = w)
          makePlot(module)
          dev.off()
          
          annRes = annotateSvgPlot(infile = svgtemp, outfile = nameimg, outdir = outdir,
            annotationInfo = module@svg, name = name)
          
          if(!annRes$annotateOK)
            svgwarn = "Note: the figure is static - enhancement with interactive effects failed. This is likely due to a version incompatibility of the 'SVGAnnotation' R package and the 'libcairo' system library. Please contact the maintainter to report this problem."
          
          sizes = paste(annRes$size)
          img = hwrite(c(aqm.hwriteImage(nameimg, width=sizes[1], height=sizes[2], id=paste("Fig", name, sep=":")),
                         annotationTable(arrayTable, name = name)))
        },
        png = {
          nameimg = paste(name, ".png", sep = "")
          png(file = file.path(outdir, nameimg), h = h*dpi, w = w*dpi)
          makePlot(module)
          dev.off()
          img = aqm.hwriteImage(nameimg, id=paste("Fig", name, sep=":"))
        },
        stop(sprintf("Invalid value of 'imageformat': %s", imageformat))
    ) ## switch

    ## Also make a PDF file
    namepdf = paste(name, ".pdf", sep = "")
    pdf(file = file.path(outdir, namepdf), h = h, w = w)
    makePlot(module)
    dev.off()
    link = list(namepdf, NA)

    ## TODO - use CSS instead of hard-coded fonts here
    linkpdf = hwrite("PDF file.",
      style='font-weight:bold;text-align:center;font-family:helvetica',
      border = 0, link = namepdf)
    
    hwrite(c(img, paste("Figure ", integerIndex, ": ", module@title,". ", linkpdf, sep="")),
           p,
           dim=c(2,1),
           style='font-weight:bold;text-align:center;font-family:helvetica',
           border=0)


    hwrite(paste("<br>", gsub("The figure <!-- FIG -->", paste("<b>Figure", integerIndex, "</b>"), module@legend, ignore.case = TRUE)), 
           p,
           style='text-align:justify;font-family:Lucida Grande;font-size:10pt')

    if(!identical(svgwarn, FALSE))
       hwrite(svgwarn,
           p,
           style='text-align:justify;font-family:Lucida Grande;font-size:10pt;color:blue')
    
  }

## End the report
makeEnding = function(p)
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


reportTable = function(p, arrayTable)
{

  ## TODO: add outlier detection

  s = seq_len(nrow(arrayTable))
  arrayTable = cbind(
    highlight =
      sprintf("<input type='checkbox' name='ReportObjectCheckBoxes' value='r:%s' onchange='checkboxEvent(\"r:%s\")'/>", s, s),
    arrayTable,
    stringsAsFactors = FALSE)

  hwrite("<hr>", p)
  hwrite(arrayTable, p, 
         row.bgcolor = rep(list("#ffffff", c("#d0d0ff", "#e0e0f0")), ceiling(nrow(arrayTable)/2)),
         table.style = "margin-left:auto;font-size:normal",
         row.style = list("font-weight:bold"))
  
}


##--------------------------------------------------
##   write the report
##--------------------------------------------------

aqm.writereport = function(modules, arrayTable, reporttitle, outdir)
  {
    ## To avoid dealing with this pathologic special case downstream in the HTML
    if(nrow(arrayTable)==0)
      stop("'arrayTable' must not be empty.")

    ## For all report modules, extract the 'svg' slot, then subset only those that are defined.
    svgdata = lapply(modules, slot, "svg")
    svgdata = svgdata[ !is.na(sapply(svgdata, slot, "name")) ]

    ## Extract strokeopacity and strokewidth from the list 'svgdata' and format for Javascript
    ##   (some of this should perhaps better be done by RJSONIO)
    formatStrokeParameters = function(w)
      paste("[", apply(sapply(svgdata, slot, w), 2, paste, collapse=", "), "]", collapse=", ")

    ## Add rownames and numeric indices to 'arrayTable'
    rown = row.names(arrayTable)
    arrayTable = if(is.numeric(rown))  ## check whether or not the row.names are 'automatic'
      cbind(row = paste(rown), arrayTable, stringsAsFactors = FALSE) else
      cbind(row = paste(seq(along=rown)), " " = paste(rown), arrayTable, stringsAsFactors = FALSE) 
    rownames(arrayTable) = NULL
             
    ## Could also use RJSONIO here
    pDataJS = sapply(arrayTable, function(x) paste("'", x, "'", sep=""))
    pDataJS = paste("[", apply(pDataJS, 1, paste, collapse=", "), "]", sep="", collapse=", ")

    ## Open and set up the HTML page
    ## Inject report-specific variables into the JavaScript
    p = makeTitle(
      reporttitle = reporttitle,
      outdir = outdir,
      params = c(          ## TODO: this should be the outliers
        HIGHLIGHTINITIAL = paste(rep(c("false", "true"), ceiling(nrow(arrayTable)/2))[seq_len(nrow(arrayTable))], collapse=", "),
        ARRAYMETADATA    = pDataJS,
        SVGOBJECTNAMES   = paste("'", names(svgdata), "'", sep="", collapse=", "),
        IDFUNS           = paste(sapply(svgdata, slot, "getPlotObjIdFromReportObjId"), collapse=", "),
        STROKEOPACITY    = formatStrokeParameters("strokeopacity"),
        STROKEWIDTH      = formatStrokeParameters("strokewidth")))

    makeIndex(p = p, modules = modules)
    reportTable(p = p, arrayTable = arrayTable)
    
    lasttype = "Something Else"
    sec = 1

    for(i in seq(along = modules))
      {
        if(modules[[i]]@section != lasttype)
          {
            makeSection(p = p, sectionNumber = sec, module = modules[[i]])
            sec = sec+1
          }
        reportModule(p = p, module = modules[[i]], integerIndex = i,
                     name = names(modules)[i], arrayTable = arrayTable, outdir=outdir)
        lasttype = modules[[i]]@section
      }
     
    makeEnding(p)
                    
    ## TODO is there a more elegant way to do this?
    invisible(list(modules=modules, arrayTable=arrayTable, reporttitle=reporttitle, outdir=outdir))
  }
