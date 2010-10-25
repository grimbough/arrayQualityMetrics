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
      if(force || length(outdirContents)==0)
        message(sprintf("The report will be written into directory '%s'. ", outdir))
        setwd(outdir)
       } else {
        dir.create(outdir, recursive=TRUE)
        message(sprintf("The directory '%s' has been created.", outdir))
        setwd(outdir)
      }
  }


## Produce the plots 
aqm.plot = function(x) {
  if (is(x@plot, "list"))
    print(x@plot) else
  if (is(x@plot, "function"))
    do.call(x@plot, args = list())
  else stop(sprintf("Invalid 'x@plot' of class '%s'.", paste(class(x@plot), collapse=", ")))
}

## Create the title
aqm.make.title = function(reporttitle, outdir)
  {
    if(!(is.character(reporttitle)&&(length(reporttitle)==1)))
      stop("'reporttitle' must be a character of length 1.")
    p = openPage(file.path(outdir, 'QMreport.html'))
    hwrite("<hr>", p)
    hwrite(reporttitle, p, heading=1, style='text-align:center;font-family:helvetica,arial,sans-serif')
    hwrite("<hr>", p)
    return(p)
  }

## Create a new section
aqm.make.section = function(p, s, qm)
  {
    hwrite("<hr>", p)
    sec = paste("<a name= 'S",s,"'>Section ", s, ": ", qm@section,"</a>", sep = "")
    hwrite(sec, p, heading=2, style='font-family:helvetica,arial,sans-serif')
  }

##To create the index
aqm.make.index = function(obj, p)
{
  s = 1
  hwrite("<hr>", p)
  hwrite("Index",p, heading=2, style='font-family:helvetica,arial,sans-serif')
  hwrite("Note: bitmap figures shown in this HTML report are also linked to PDF files, to provide better figure quality, or to provide multi-page files in cases where in  this HTML report only the first page is presented.", p, style='font-family:helvetica;font-size:11pt;color:#808080;font-style:italic;')

  hwrite("<UL>", p)
  lasttype = "FAKE"
  for(i in seq_len(length(obj)))
    {
      if(obj[[i]]$section != lasttype)
        {
          if(s != 1)
            hwrite("</UL>", p)

          hwrite(paste("<br><LI>Section ", s,": ", obj[[i]]$section,"<UL>",sep=""), p, link=paste("#S",s,sep=""), style= 'font-weight:bold;font-family:helvetica;font-size:12pt')
          s = s+1
        }
      hwrite(paste("<LI>",obj[[i]]$title,sep=""), p, style= 'font-weight:normal;font-family:helvetica;font-size:11pt')
      lasttype = obj[[i]]$section
    }
  hwrite("</UL></UL>", p)
}


##To create a part of report with figures and legend
aqm.report.qm = function(p, qm, f, name, outdir)
  {
    stopifnot(is(qm, "aqmReportModule"))

    h = qm@shape$h
    w = qm@shape$w
    dpi = 72

    imageformat =  if(length(qm@svg)>0) "svg" else "png"

    ## The two cases, png and svg, need to be treated differently 
    switch(imageformat,
           svg = {
             nameimg = paste(name, ".svg", sep = "")
             svgtemp = tempfile()
             svg(file = svgtemp, h = h, w = w)
             aqm.plot(qm)
             dev.off()
             size = annotateSvgMatplot(svgtemp, file.path(outdir, nameimg), annotationInfo=qm@svg)
             img = aqm.hwriteImage(nameimg, width=paste(size[1]), height=paste(size[2]))
           },
           png = {
             nameimg = file.path(outdir, paste(name, ".png", sep = ""))
             png(file = nameimg, h = h*dpi, w = w*dpi)
             aqm.plot(qm)
             dev.off()
             img = aqm.hwriteImage(nameimg)
           },
           stop(sprintf("Invalid value of 'imageformat': %s", imageformat))
           )
    namepdf = paste(name, ".pdf", sep = "")
    pdf(file = namepdf, h = h, w = w)
    aqm.plot(qm)
    dev.off()
    link = list(namepdf, NA)
      
    linkpdf = hwrite("PDF file.",
      style='font-weight:bold;text-align:center;font-family:helvetica',
      border=0, link=namepdf)
    
    hwrite(c(img, paste("Figure ", f, ": ", qm@title,". ", linkpdf, sep="")),
           p,
           dim=c(2,1),
           style='font-weight:bold;text-align:center;font-family:helvetica',
           border=0, link=link)
    
    hwrite(paste("<br>", qm@legend),
           p,
           style='text-align:justify;font-family:Lucida Grande;font-size:10pt')
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
scores = function(expressionset, obj)
  {
    titles = c(
      "aqmobj.ma"      = "MA plots",
      "aqmobj.spatial" = "Spatial",
      "aqmobj.box"     = "Distribution",
      "aqmobj.heat"    = "Heatmap",
      "aqmobj.rle"     = "RLE",
      "aqmobj.nuse"    = "NUSE")

    classes = sapply(obj, class)
    mt = match(classes, names(titles))

    sN = sampleNames(expressionset)

    df = data.frame(
      "Array #"    = seq_len(length(sN)),
      "Array Name" = sN,
      stringsAsFactors = FALSE)
  
    if(!is.null(protocolData(expressionset)$ScanDate))
      df$"Scan Dates" = protocolData(expressionset)$ScanDate
    
    for(i in which(!is.na(mt)))      
      df[[titles[mt[i]]]] =
        ifelse(seq(along=sN) %in% obj[[i]]$outliers, "*", "")

    return(df)
  }


aqm.writereport = function(name, obj, outdir)
  {
    ## TODO: do we need this? Should it really be done here, or next to the processing of the 'usesvg' argument
    ##    in the 'arrayQualityMetrics' function?
    if (Sys.info()["sysname"] == "Darwin")
      options(bitmapType = "cairo")

    sec = 1
    p = aqm.make.title(name, outdir)
  
    sc = scores(obj)
    col = rep(c("#d0d0ff","#e0e0f0"), (ceiling((nrow(sc)+1)/2))) #

    if(nrow(sc) < length(col))  
      col = col[seq_len((length(col)-(abs((nrow(sc)+1) - length(col)))))]

    boldset = 'font-weight:bold;text-align:center;font-family:Lucida Grande;font-size:10pt'
    normset = 'text-align:center;font-family:Lucida Grande;font-size:10pt'

    hwrite("Summary",p, heading=2, style='font-family:helvetica,arial,sans-serif')
    hwrite(sc, p, border=0, bgcolor = matrix(col, ncol=ncol(sc), nrow=(nrow(sc)+1)), cellpadding = 2, cellspacing = 5, style=matrix(c(rep(boldset, ncol(sc)), rep(c(boldset, rep(normset, ncol(sc)-1)), nrow(sc))), ncol=ncol(sc), nrow=(nrow(sc)+1), byrow=TRUE))
          
    hwrite("*outlier array",style=normset,p)


    aqm.make.index(obj, p)

    lasttype = "Something Else"
    
    for(i in seq(along = obj))
      {
        if(obj[[i]]$section != lasttype)
          {
            aqm.make.section(p, s = sec, qm = obj[[i]])
            sec = sec+1
          }
        aqm.report.qm(p, qm=obj[[i]], f=i, name=names(obj)[i], outdir=outdir)
        lasttype = obj[[i]]$section
      }
     
    aqm.make.ending(p)
  }
