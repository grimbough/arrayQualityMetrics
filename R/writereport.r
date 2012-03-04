##---------------------------------------------------------
## Create the output directory
##---------------------------------------------------------
dircreation = function(outdir = getwd(), force = FALSE)
  {
    if(file.exists(outdir)){
      if(!file.info(outdir)$isdir)
        stop(sprintf("'%s' must be a directory.", outdir))

      outdirContents = dir(outdir, all.files = TRUE)
      outdirContents = setdiff(outdirContents, c(".", ".."))

      if(!force && length(outdirContents)>0)
        stop(sprintf("The directory '%s' already exists and is not empty. Please remove the directory before calling 'arrayQualityMetrics' or consider using 'force=TRUE'.", outdir))
      message(sprintf("The report will be written into directory '%s'. ", outdir))
    } else {
      dir.create(outdir, recursive=TRUE)
      message(sprintf("The directory '%s' has been created.", outdir))
    }
  }


##---------------------------------------------------------
## Produce a plot
##---------------------------------------------------------
makePlot = function(x) {
  if (is(x@plot, "trellis"))
    print(x@plot) else
  if (is(x@plot, "function"))
    do.call(x@plot, args = list())
  else stop(sprintf("Invalid 'x@plot' of class '%s'.", paste(class(x@plot), collapse=", ")))
}

##---------------------------------------------------------
## Create the title
##---------------------------------------------------------
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

    p = openPage(filename = file.path(outdir, 'index.html'),
      link.javascript = filenames[2],
      link.css        = filenames[1],
      body.attributes = list("onload" = "reportinit()"),
      title           = reporttitle)

    hwrite("<hr>", page = p)
    hwrite(reporttitle, page = p, heading=1)
    hwrite("<hr>", page = p)
    return(p)
  }


##---------------------------------------------------------
## Create a new section
##---------------------------------------------------------
makeSection = function(p, sectionNumber, module)
{
  hwrite("<hr>", page = p)
  sec = paste("<a name= 'S", sectionNumber,"'>Section ", sectionNumber, ": ",
    module@section,"</a>", sep = "")
  hwrite(sec, page = p, heading=2)
}

##---------------------------------------------------------
## Create the index
##---------------------------------------------------------
makeIndex = function(p, modules)
{
  currentSectionNumber = 1
  currentSectionName   = "Something else"

  hwrite("<UL>", page = p)
  for(i in seq_len(length(modules)))
    {
      if(modules[[i]]@section != currentSectionName)
        {
          if(currentSectionNumber != 1) ## end the previous section
            hwrite("</UL>", page = p)

          hwrite(paste("<br><li class='tocsection'>Section ", currentSectionNumber,": ",
                       modules[[i]]@section, "</li><UL>", sep=""), page = p,
                 link = paste("#S", currentSectionNumber, sep=""))
          currentSectionNumber = currentSectionNumber+1
        }
      hwrite(paste("<li class='tocmodule'>", modules[[i]]@title, "</li>", sep=""), page = p)
      currentSectionName = modules[[i]]@section
    }
  hwrite("</UL></UL>", page = p)

  browserCompatibilityNote(p)
}

##---------------------------------------------------------
## Create a module  of the report with figures and legend
##---------------------------------------------------------
reportModule = function(p, module, currentIndex, arrayTable, outdir)
{
    stopifnot(is(module, "aqmReportModule"))
    validObject(module, test=FALSE)

    svgwarn = FALSE

    stopifnot(!is.na(module@title))
    name = module@id

    stopifnot(!any(is.na(module@size)))
    h = module@size["h"]
    w = module@size["w"]

    stopifnot(length(currentIndex)==1, is.numeric(currentIndex), !is.na(currentIndex))
    dpi = arrayQualityMetricsGlobalParameters$dpi

    if(is.na(module@svg@numPlotObjects))
      {
        ## no svg - use png
        nameimg = paste(name, ".png", sep = "")
        png(filename = file.path(outdir, nameimg), height = h*dpi, width = w*dpi)
        makePlot(module)
        dev.off()
        img = hmakeTag("img", src = nameimg, border = 0,
                       alt = nameimg, id = paste("Fig", name, sep=":"))

      } else {
        ## svg
        nameimg = paste(name, ".svg", sep = "")
        svgtemp = paste(tempfile(), ".svg", sep = "")
        Cairo(file = svgtemp, type = "svg", height = h, width = w, units = "in", dpi = dpi)
        makePlot(module)
        dev.off()

        annRes = annotateSvgPlot(infile = svgtemp,
                                 outfile = nameimg,
                                 outdir = outdir,
                                 annotationInfo = module@svg,
                                 name = name)

        if(!annRes$annotateOK)
          svgwarn = paste("Note: the figure is static - enhancement with interactive effects failed.",
            "This is either due to a version incompatibility of the 'SVGAnnotation' R package and your",
            "version of 'Cairo' or 'libcairo', or due to plot misformating. Please consult the Bioconductor mailing list, or",
            "contact the maintainer of 'arrayQualityMetrics' with a reproducible example in order to fix this problem.")

        sizes = paste(round(annRes$size))
        img = hwrite( c( paste(readLines(file.path(outdir, nameimg)), collapse="\n"),
                         annotationTable(arrayTable, name = name) ))

        ## TO DO:
        ##  hwrite(c(aqm.hwriteImage(nameimg, width=sizes[1], height=sizes[2], id=paste("Fig", name, sep=":")),
        ##               annotationTable(arrayTable, name = name)))
      } ## if

    ## Also make a PDF file
    namepdf = paste(name, ".pdf", sep = "")
    pdf(file = file.path(outdir, namepdf), height = h, width = w)
    makePlot(module)
    dev.off()

    # write the HTML
    hwrite("\n\n", page = p)
    hwrite(toggleStart(name, display=module@defaultdisplay, text = sprintf("Figure %d: %s.", currentIndex, module@title)), page = p)

    hwrite(img, page = p)
    hwrite("<br>\n", page = p)

    hwrite(gsub("The figure <!-- FIG -->",
           paste("<b>Figure ", currentIndex, "</b>", if(!is.na(namepdf)) hwrite(" (PDF file)", link = namepdf), sep=""),
                 module@legend, ignore.case = TRUE), page = p)
    hwrite("<br><br><br>\n", page = p)

    if(!identical(svgwarn, FALSE))
       hwrite(svgwarn, page = p)

    hwrite(toggleEnd(), page = p)

    ## recursion, for the barplot with the outliers
    if(!identical(NA_character_, module@outliers@description)) {
      currentIndex = currentIndex + 1
      reportModule(p, aqm.outliers(module), currentIndex, arrayTable, outdir)
    }

    return(currentIndex + 1)
  }

##----------------------------------------------------------
## End the report
##----------------------------------------------------------
makeEnding = function(p)
  {
    z = sessionInfo("arrayQualityMetrics")
    version = z$otherPkgs[[1]]$Version
    rversion = sessionInfo()$R.version$version.string
    session = paste("This report has been created with arrayQualityMetrics ",
      version, " under ", rversion, ".", sep="")
    hwrite("<hr>", page = p)
    hwrite(session, page = p, style ='font-size:9pt')
    hwrite("<hr>", page = p)
    closePage(page = p)
  }

##----------------------------------------------------------
## Browser compatibility
##----------------------------------------------------------
browserCompatibilityNote = function(p)
  {
    txt = paste("<h3>Browser compatibility</h3>\n",
      "This report uses recent features of HTML 5. Functionality has been tested on these browsers: ",
      "Firefox 10, Chrome 17, Safari 5.1.2\n", sep="")
    hwrite("<hr>\n", page = p)
    hwrite(txt, page = p)
  }

##----------------------------------------------------------
##   write the HTML for 'arrayTable', including checkboxes
##----------------------------------------------------------
reportTable = function(p, arrayTable, tableLegend)
{
  s = seq_len(nrow(arrayTable))
  arrayTable = cbind(
    " " = sprintf("<input type='checkbox' name='ReportObjectCheckBoxes' value='' onchange='checkboxEvent(%d)'>", s),
    arrayTable,
    stringsAsFactors = FALSE)

  hwrite("<hr>", page = p)

  disp = ifelse(nrow(arrayTable)<=arrayQualityMetricsGlobalParameters$maxNumberOfArraysForShowingArrayMetadataByDefault , "block", "none")
  hwrite(toggleStart("arraymetadata", disp, "Array metadata and outlier detection overview"), page = p)
  hwrite(arrayTable, page = p,
         row.bgcolor = rep(list("#ffffff", c("#d0d0ff", "#e0e0f0")), ceiling(nrow(arrayTable)/2)),
         table.style = "margin-left:auto;text-align:right;",
         row.style = list("font-weight:bold"))

  hwrite(paste("<br>", tableLegend, "<br>", toggleEnd(), sep=""), page = p)
}

##------------------------------------------------------------------
## Display toggle start and end
##----------------------------------------------------------
toggleStart = function(name, display, text)
  paste(
    sprintf("<a name=\"%s\" id=\"%s-h\" href=\"javascript: toggle('%s')\" style=\"text-decoration:none;font-weight:bold;font-size:larger\">%s %s</a><br>\n",
                      name,      name,                            name,                                  c("block"="-", "none"="+")[display], text),
    sprintf("<div id=\"%s-b\" style=\"display:%s\">\n",
                      name,                 display), sep="\n")

toggleEnd = function()
  "</div>"

##------------------------------------------------------------------
## Create JSON representation of R character vectors and matrices
## Names and dimnames attributes are stripped.
##----------------------------------------------------------
toJSON_fromchar = function(x)
  paste("[", paste(x, collapse=", "), "]")

toJSON_fromvector = function(x)
    toJSON_fromchar(paste('"', as.character(x), '"', sep=''))

toJSON_frommatrix = function(x)
  {
    stopifnot(length(dim(x))==2)
    toJSON_fromchar(apply(x, 1, toJSON_fromvector))
  }

##--------------------------------------------------
##   write the report
##--------------------------------------------------
aqm.writereport = function(modules, arrayTable, reporttitle, outdir)
{
  numReportObjs = nrow(arrayTable)
  reportObjs    = seq_len(numReportObjs)

  ## To avoid dealing with this pathologic special case downstream in the HTML
  if(numReportObjs==0)
    stop("'arrayTable' must not be empty.")

  ## construct short, unique IDs
  ids = sapply(modules, slot, "id")
  stopifnot(!any(is.na(ids)), !any(duplicated(modules)))

  ## For all report modules, extract the 'svg' slot, then subset only those that are defined.
  svgdata = lapply(modules, slot, "svg")
  names(svgdata) = ids
  hassvg  = !is.na(sapply(svgdata, slot, "numPlotObjects"))
  svgdata = svgdata[ hassvg]

  ## Determine which subset of the modules have computed outliers ('wh'). For each, define
  ## a corresponding column in the logical matrix 'outliers'.
  ## Further below, a textual representation of 'outliers', ifelse(outliers, "x", "") is added to arrayTable,
  ## and the row-wise OR (more precisely: 'apply(outliers, 1, any)') is used to determine
  ## which arrays to highlight initially.
  wh = which(sapply(modules, function(x) length(x@outliers@statistic)>0))

  outlierMethodTitles = sapply(modules, slot, "title")[wh]
  outlierMethodLinks  = paste("<a href=\"#", ids[wh], "\">", sep="")

  outlierExplanations = paste(
    "The columns named *1, *2, ... indicate the calls from the different outlier detection methods:<OL>",
    paste(sprintf("<LI> outlier detection by %s%s</a></LI>",
                       outlierMethodLinks, outlierMethodTitles), collapse = ""),
    "</OL>The outlier detection criteria are explained below in the respective sections. Arrays that were called outliers ",
    "by at least one criterion are marked by checkbox selection in this table, and are ",
    "indicated by highlighted lines or points in some of the plots below. ",
    "By clicking the checkboxes in the table, or on the corresponding points/lines in the plots, you can modify the selection. ",
    "To reset the selection, reload the HTML page in your browser.", sep="")

  outliers = matrix(NA, nrow = numReportObjs,
                        ncol = length(wh),
                        dimnames = list(NULL, sprintf("%s*%d</a>", outlierMethodLinks, seq(along=wh))))

  for(j in seq(along = wh))
    {
      o = modules[[wh[j]]]@outliers@which
      stopifnot(!any(is.na(o)), all( (o>=1) & (o<=numReportObjs)))
      outliers[,j] = reportObjs %in% o
    }


  ## Add numeric indices, rownames and outlier annotation to 'arrayTable'
  ## Make two versions of it:
  ## - 'big', includes outlier status, is shown at the top of the report
  ## - 'compact' , without outlier status, is shown next to the interactive plots
  rowchar = as.character(row.names(arrayTable))
  rownum  = paste(reportObjs)
  left = if(!identical(rowchar, rownum))
    data.frame(array = rownum, sampleNames = rowchar, stringsAsFactors = FALSE) else
    data.frame(array = rownum, stringsAsFactors = FALSE)

  arrayTableBig     = cbind(left, ifelse(outliers, "x", ""), arrayTable, stringsAsFactors = FALSE)
  arrayTableCompact = cbind(left, arrayTable, stringsAsFactors = FALSE)
  rownames(arrayTableBig) = rownames(arrayTableCompact) = NULL

  ## Open and set up the HTML page
  p = makeTitle(
    reporttitle = reporttitle,
    outdir = outdir,
    ## Inject report-specific variables into the JavaScript
    params = c(
      HIGHLIGHTINITIAL = toJSON_fromchar(ifelse(apply(outliers, 1, any), "true", "false")),
      ARRAYMETADATA    = toJSON_frommatrix(arrayTableCompact),
      SVGOBJECTNAMES   = toJSON_fromvector(names(svgdata)),
      REPORTOBJSTYLES  = paste(".aqm", reportObjs, " { }", sep="", collapse = "\n")
    ))

  makeIndex(p = p, modules = modules)
  reportTable(p = p, arrayTable = arrayTableBig,
              tableLegend = outlierExplanations)

  currentSectionName = "Something Else"
  currentIndex = currentSection = 1

  for(i in seq(along = modules))
    {
      if(modules[[i]]@section != currentSectionName)
        {
          makeSection(p = p, sectionNumber = currentSection, module = modules[[i]])
          currentSection = currentSection+1
        }
      currentIndex = reportModule(
        p = p,
        module = modules[[i]],
        currentIndex = currentIndex,
        arrayTable = arrayTableCompact,
        outdir=outdir)
      currentSectionName = modules[[i]]@section
    }

  makeEnding(p)
  invisible(list(modules=modules, arrayTable=arrayTableBig, reporttitle=reporttitle, outdir=outdir))
}
