setClassUnion("aqmTrellis", c("aqmobj.ma", "aqmobj.heat", "aqmobj.qcs", "aqmobj.probesmap", "aqmobj.spatial", "aqmobj.spatialbg"))

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
  if(class(obj$plot) == "trellis")
    print(obj$plot)
  if(class(obj$plot) == "list")
    {
      print(obj$plot$boxred, split = c(1,1,3,1), newpage = FALSE)
      print(obj$plot$boxgreen, split = c(2,1,3,1), newpage = FALSE)
      print(obj$plot$boxblue, split = c(3,1,3,1), newpage = FALSE)
    }
})

setMethod("aqm.plot",signature(obj = "aqmobj.dens"), function(obj){
  if(class(obj$plot) == "trellis")
    print(obj$plot)
  if(class(obj$plot) == "list")
  {
    print(obj$plot$den1, split = c(1,1,3,1), newpage = FALSE)
    print(obj$plot$den2, split = c(2,1,3,1), newpage = FALSE)
    print(obj$plot$den3, split = c(3,1,3,1), newpage = FALSE)
  }
})

setMethod("aqm.plot",signature(obj = "aqmobj.rle"), function(obj){
  box.rectangle <- trellis.par.get("box.rectangle")
  box.rectangle$col = "black"
  trellis.par.set("box.rectangle",box.rectangle)
  box.umbrella <- trellis.par.get("box.umbrella")
  box.umbrella$col = "transparent"
  trellis.par.set("box.umbrella",box.umbrella)
  bwplot(obj$plot$stats ~ as.vector(col(obj$plot$stats)), pch = "|", col = "black", do.out = FALSE, fill = "#1F78B4", horizontal = FALSE, box.ratio = 2, ylim = c(-1,1), ylab="RLE", xlab = "")
})

setMethod("aqm.plot",signature(obj = "aqmobj.nuse"), function(obj){
  box.rectangle <- trellis.par.get("box.rectangle")
  box.rectangle$col = "black"
  trellis.par.set("box.rectangle",box.rectangle)
  box.umbrella <- trellis.par.get("box.umbrella")
  box.umbrella$col = "black"
  trellis.par.set("box.umbrella",box.umbrella)
  bwplot(obj$plot$stats ~ as.vector(col(obj$plot$stats)), pch = "|", col = "black", do.out = FALSE, fill = "#1F78B4", horizontal = FALSE, box.ratio = 2, ylim = c(0.9,1.5), ylab="NUSE", xlab = "")  
})

setMethod("aqm.plot",signature(obj = "aqmobj.rnadeg"), function(obj){
  acol = sample(brewer.pal(8, "Dark2"), length(obj$plot$sample.names), replace = (8<length(obj$plot$sample.names)))
  plotAffyRNAdeg(obj$plot, lwd = 2, col =acol)})

setMethod("aqm.plot",signature(obj = "aqmobj.msd"), function(obj){
  meanSdPlot(obj$plot)})
  
##To create the title
aqm.make.title = function(arg)
  {
    argum = arg$expressionset
    p = openPage('QMreport.html')
    hwrite("<hr>", p)
    title = paste(deparse(substitute(argum)), " quality metrics report", sep="")
    hwrite(title, p, heading=1, style='text-align:center;font-family:helvetica,arial,sans-serif')
    hwrite("<hr>", p)
    return(p)
  }

##To create a new section
aqm.make.section = function(p, s, qm)
  {
    hwrite("<hr>", p)
    sec = paste("<a name= 'S",s,"'>Section ", s, ": ", qm$type,"</a>", sep = "")
    hwrite(sec, p, heading=2, style='font-family:helvetica,arial,sans-serif')
  }


##To create the index
aqm.make.index = function(obj, p)
{
  hwrite("<hr>", p)
  hwrite("Index",p, heading=2, style='font-family:helvetica,arial,sans-serif')

  hwrite("<UL>", p)
  s = 1
  lasttype = "FAKE"
  for(i in seq_len(length(obj)))
    {
      if(obj[[i]]$type != lasttype)
        {
          if(s != 1)
            hwrite("</UL>", p)

          hwrite(paste("<br><LI>Section ", s,": ", obj[[i]]$type,"<UL>",sep=""), p, link=paste("#S",s,sep=""), style= 'font-weight:bold;font-family:helvetica;font-size:12pt')
          s = s+1
        }
      hwrite(paste("<LI>",obj[[i]]$title,sep=""), p, style= 'font-weight:normal;font-family:helvetica;font-size:11pt')
      lasttype = obj[[i]]$type
    }
  hwrite("</UL></UL>", p)
}


##To create a part of report with figures and legend
aqm.report.qm = function(p, qm, f, name)
  {
    if(qm$shape == "rect")
      {
        h = 5
        w = 15
      }
    if(qm$shape == "square")
      {
        h = 7
        w = 7
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
    hwrite(c(qm$title, "", img,  paste("Figure",f)), p, dim=c(2,2), style='font-weight:bold;text-align:center;font-family:helvetica', border=0, link=list(NA,NA,namepdf,namepdf))
    hwrite(qm$legend, p, style='text-align:justify;font-family:Lucida Grande;font-size:10pt')
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
    if(inherits(expressionset, "BeadLevelList"))
      sN = arrayNames(expressionset)
    if(inherits(expressionset, "NChannelSet"))
      sN = sampleNames(expressionset)$R
    if(inherits(expressionset, "ExpressionSet") || inherits(expressionset, "AffyBatch"))
      sN = sampleNames(expressionset)

    nscores = 0
    namesm = c()

    if(length(obj$maplot) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "MA plots")
      }
    if(length(obj$spatial) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "Spatial distribution")
      }
    if(length(obj$boxplot) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "Boxplots")
      }
    if(length(obj$heatmap) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "Heatmap")
      }
    if(length(obj$rle) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "RLE")
      }
    if(length(obj$nuse) != 0)
      {
        nscores = nscores +1
        namesm = c(namesm, "NUSE")
      }
  
    m = matrix("",nrow = length(sN), ncol = nscores+2)
    colnames(m) = c("Array #", "Array Name", namesm)
    m[,1] = seq_len(length(sN))
    m[,2] = sN
    
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


