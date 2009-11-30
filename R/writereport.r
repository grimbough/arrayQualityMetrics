setClassUnion("aqmTrellis", c("aqmobj.ma", "aqmobj.heat", "aqmobj.pca", "aqmobj.probesmap", "aqmobj.spatial", "aqmobj.spatialbg", "aqmobj.dens"))

##Creation of the outdir
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
        message(sprintf("The report will be written in directory '%s'. ", outdir))
        setwd(outdir)
       } else {
        dir.create(outdir, recursive=TRUE)
        message(sprintf("The directory '%s' has been created.", outdir))
        setwd(outdir)
      }
  }

## Write report functions
setGeneric("aqm.plot",
           function(obj)
           standardGeneric("aqm.plot"))

##To produce the plots
setMethod("aqm.plot",signature(obj = "aqmTrellis"), function(obj){
  print(obj$plot)})

setMethod("aqm.plot",signature(obj = "aqmobj.pmmm"), function(obj){
plot(obj$plot$MM, col = "grey", xlab = "log(Intensity)", main="")
lines(obj$plot$PM, col = "blue")
legend("topright", c("PM","MM"),lty=1,lwd=2,col= c("blue","grey"), bty = "n")
})

setMethod("aqm.plot",signature(obj = "aqmobj.box"), function(obj){
  box.rectangle <- trellis.par.get("box.rectangle")
  box.rectangle$col = "black"
  trellis.par.set("box.rectangle",box.rectangle)
  box.umbrella <- trellis.par.get("box.umbrella")
  box.umbrella$col = "black"
  trellis.par.set("box.umbrella",box.umbrella)
  print(obj$plot)
})

setMethod("aqm.plot",signature(obj = "aqmobj.rle"), function(obj){
  box.rectangle <- trellis.par.get("box.rectangle")
  box.rectangle$col = "black"
  trellis.par.set("box.rectangle",box.rectangle)
  box.umbrella <- trellis.par.get("box.umbrella")
  box.umbrella$col = "transparent"
  trellis.par.set("box.umbrella",box.umbrella)
  print(bwplot(obj$plot$stats ~ as.vector(col(obj$plot$stats)), pch = "|", col = "black", do.out = FALSE, fill = "#1F78B4", horizontal = FALSE, box.ratio = 2, ylim = c(-1,1), ylab="RLE", xlab = ""))
})

setMethod("aqm.plot",signature(obj = "aqmobj.nuse"), function(obj){
  box.rectangle <- trellis.par.get("box.rectangle")
  box.rectangle$col = "black"
  trellis.par.set("box.rectangle",box.rectangle)
  box.umbrella <- trellis.par.get("box.umbrella")
  box.umbrella$col = "black"
  trellis.par.set("box.umbrella",box.umbrella)
  print(bwplot(obj$plot$stats ~ as.vector(col(obj$plot$stats)), pch = "|", col = "black", do.out = FALSE, fill = "#1F78B4", horizontal = FALSE, box.ratio = 2, ylim = c(0.9,1.5), ylab="NUSE", xlab = ""))  
})

setMethod("aqm.plot",signature(obj = "aqmobj.rnadeg"), function(obj){
  acol = sample(brewer.pal(8, "Dark2"), length(obj$plot$sample.names), replace = (8<length(obj$plot$sample.names)))
  print(plotAffyRNAdeg(obj$plot, lwd = 2, col =acol))
  legend("topright",lty=1,lwd=2,col=acol,legend = obj$plot$sample.names)


})

setMethod("aqm.plot",signature(obj = "aqmobj.qcs"), function(obj){
  plot.qc.stats(obj$plot)})
  
setMethod("aqm.plot",signature(obj = "aqmobj.msd"), function(obj){
  meanSdPlot(obj$plot)})
  
##To create the title
aqm.make.title = function(name)
  {
    p = openPage('QMreport.html')
    hwrite("<hr>", p)
    title = paste(name, " quality metrics report", sep="")
    hwrite(title, p, heading=1, style='text-align:center;font-family:helvetica,arial,sans-serif')
    hwrite("<hr>", p)
    return(p)
  }

##To create a new section
aqm.make.section = function(p, s, qm)
  {
    hwrite("<hr>", p)
    sec = paste("<a name= 'S",s,"'>Section ", s, ": ", qm$section,"</a>", sep = "")
    hwrite(sec, p, heading=2, style='font-family:helvetica,arial,sans-serif')
  }


##To create the index
aqm.make.index = function(obj, p)
{
  s = 1
  hwrite("<hr>", p)
  hwrite("Index",p, heading=2, style='font-family:helvetica,arial,sans-serif')
  hwrite("PLEASE NOTE:<br>All figures below are links to PDF files: these contain images for every array in the report.  The PDF files may be several pages long, this HTML report presents only the first page.", p, style='font-weight:bold;font-family:helvetica;font-size:11pt;color:#FF0000')

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
aqm.report.qm = function(p, qm, f, name)
  {
    if(qm$shape == "rect")
      {
        h = 6
        w = 10
      }
    if(qm$shape == "square")
      {
        h = 6
        w = 6
      }
    
    dpi = 72

    namepdf = paste(name, ".pdf", sep = "")
    namepng = paste(name, ".png", sep = "")

    png(file = namepng, h = h*dpi, w = w*dpi)
    if(class(qm) == "aqmobj.ma"  || ((class(qm) == "aqmobj.spatial" || class(qm) == "aqmobj.spatialbg") && (class(qm$plot) == "list" && class(qm$plot[[1]]) == "trellis")))
      print(qm$plot[[1]]) else aqm.plot(qm)
    dev.off()

    pdf(file = namepdf, h = h, w = w)
    aqm.plot(qm)
    dev.off()   
  
    img = hwriteImage(namepng)
    hwrite(c(img,  paste("Figure ",f,": ",qm$title, sep="")), p, dim=c(2,1), style='font-weight:bold;text-align:center;font-family:helvetica', border=0, link=list(namepdf,NA))
    hwrite(paste("<br>",qm$legend), p, style='text-align:justify;font-family:Lucida Grande;font-size:10pt')
  }

##To end the report
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

##Score table formatting
scores = function(expressionset, obj)
  {
    sN = sampleNames(expressionset)
    nscores = 0
    namesm = c()

    classes = sapply(seq(along = obj), function(i) class(obj[[i]]))

    if(length(grep("aqmobj.ma",classes)) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "MA plots")
      }
    if(length(grep("aqmobj.spatial",classes)) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "Spatial distribution")
      }
    if(length(grep("aqmobj.box",classes)) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "Boxplots")
      }
    if(length(grep("aqmobj.heat",classes)) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "Heatmap")
      }
    if(length(grep("aqmobj.rle",classes)) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "RLE")
      }
    if(length(grep("aqmobj.nuse",classes)) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "NUSE")
      }
  
    if(is.null(protocolData(expressionset)$ScanDate))
    {
    m = matrix("",nrow = length(sN), ncol = nscores+2)
    colnames(m) = c("Array #", "Array Name", namesm)
    m[,1] = seq_len(length(sN))
    m[,2] = sN    
    } else {
    m = matrix("",nrow = length(sN), ncol = nscores+3)
    colnames(m) = c("Array #", "Array Name", "Scan Dates", namesm)
    m[,1] = seq_len(length(sN))
    m[,2] = sN
    m[,3] = protocolData(expressionset)$ScanDate    
    }
    
    if(length(obj$maplot$outliers) != 0)
      m[obj$maplot$outliers,"MA plots"] = "*"    
    if(length(obj$spatial$outliers) != 0)
      m[unlist(obj$spatial$outliers),"Spatial distribution"] = "*"    
    if(length(unlist(obj$boxplot$outliers)) != 0)
      m[unlist(obj$boxplot$outliers),"Boxplots"] = "*"    
    if(length(obj$heatmap$outliers) != 0)
      m[obj$heatmap$outliers,"Heatmap"] = "*"    
    if(length(obj$rle$outliers) != 0)
      m[obj$rle$outliers,"RLE"] = "*"    
    if(length(unlist(obj$nuse$outliers)) != 0)
      m[unlist(obj$nuse$outliers),"NUSE"] = "*"

    return(m)
  }



aqm.writereport = function(name, expressionset, obj)
  {
    if (Sys.info()["sysname"] == "Darwin")
      options(bitmapType = "cairo")
    sec = 1
    p = aqm.make.title(name)
  
    sc = scores(expressionset, obj)
    col = rep(c("#d0d0ff","#e0e0f0"), (ceiling((nrow(sc)+1)/2))) #

    if(nrow(sc) < length(col))  
      col = col[seq_len((length(col)-(abs((nrow(sc)+1) - length(col)))))]

    boldset = 'font-weight:bold;text-align:center;font-family:Lucida Grande;font-size:10pt'
    normset = 'text-align:center;font-family:Lucida Grande;font-size:10pt'

    hwrite("Summary",p, heading=2, style='font-family:helvetica,arial,sans-serif')
    hwrite(sc, p, border=0, bgcolor = matrix(col, ncol=ncol(sc), nrow=(nrow(sc)+1)), cellpadding = 2, cellspacing = 5, style=matrix(c(rep(boldset, ncol(sc)), rep(c(boldset, rep(normset, ncol(sc)-1)), nrow(sc))), ncol=ncol(sc), nrow=(nrow(sc)+1), byrow=TRUE))
          
    hwrite("*outlier array",style=normset,p)

    if(class(obj) != "list")
      {
        aqm.make.section(p, s = sec, qm = obj)
        aqm.report.qm(p, obj, 1, "fig")
      } else {

        aqm.make.index(obj, p)

        lasttype = "FAKE"
    
        for(i in seq(along = obj))
          {
            if(obj[[i]]$section != lasttype)
              {
                aqm.make.section(p, s = sec, qm = obj[[i]])
                sec = sec+1
              }
            aqm.report.qm(p, obj[[i]], i, names(obj)[i])
            lasttype = obj[[i]]$section
          }
      }
    
    aqm.make.ending(p)
  }
