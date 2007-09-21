setGeneric("arrayQualityMetrics",
           function(expressionset,
                    outdir = getwd(),
                    force = FALSE,
                    do.logtransform = FALSE,
                    split.plots = FALSE)
           standardGeneric("arrayQualityMetrics"))

#################################################################
#################################################################
########################### FUNCTIONS ###########################
#################################################################
#################################################################

##lists
mat2list = function(x)
      lapply(seq_len(ncol(x)), function(i) x[,i])

##makePlot
makePlot = function(con, name, w, h=devDims(w)$height, fun, psz=12, isPlatePlot=FALSE, isImageScreen=FALSE, title, text, fig) {
  outf = paste(name, c("pdf", "png"), sep=".")
  nrppi = 72

  pdf(outf[1], width=w, height=h, pointsize=psz)
  if (isImageScreen)  fun(map=FALSE) else fun()
  dev.off()

  if (isPlatePlot) {
    wd = devDims(w)$pwidth
    hg = devDims(w)$pheight
  } else {
    wd = w*nrppi
    hg = h*nrppi
  }
  
  png(outf[2], width=wd, height=hg, pointsize=psz)
  res <- fun()
  dev.off()
  
  htmltext = sprintf(text, title, basename(outf[1]), basename(outf[2]), fig)  
  resret = list(res,htmltext)  
  return(resret)
}

##Creation of the outdir
dircreation = function(outdir = getwd(), force = FALSE)
  {
    if(file.exists(outdir)){
      if(!file.info(outdir)$isdir)
        stop(sprintf("'%s' must be a directory.", outdir))
      outdirContents = dir(outdir, all.files = TRUE)
      outdirContents = setdiff(outdirContents, c(".", ".."))
              
      if(length(outdirContents)>0) {
        if(!force)
          stop(sprintf("'%s' is not empty.", outdir))
        unlink(file.path(outdir, outdirContents), recursive=TRUE)
      } 
    } else {
      dir.create(outdir, recursive=TRUE)
      message(sprintf("The directory '%s' has been created.", outdir))
    }
    setwd(outdir)
  }

##second part of data preparation
prepdata = function(sN, dat, numArrays, split.plots)
  {
    ##maximum length of the experiment names to adjust margins and font size on the plot axis
    long = max(nchar(sN))
      
    ##if experiment names are to long they are replaced by figures
    if(long >= 20)
      {
        sNt = cbind(sN,seq_len(length(sN)))        
        colnames(sNt)=c("SampleName","New denomination")
        sN = seq_len(length(sN))
        long = max(nchar(sN))
      }
    if(numArrays >= 50 && long >= 4)
      {
        sNt = cbind(sN,seq_len(length(sN)))        
        colnames(sNt)=c("SampleName","New denomination")
        sN = seq_len(length(sN))
        long = max(nchar(sN))
      }

    outM = as.dist(dist2(na.omit(dat)))
      
    k = if(split.plots) split.plots else k = numArrays
    ##attribute randomly the experiments to different groups
    group = sample(rep((1:ceiling(numArrays/k)),k),numArrays)
    
    dp = if(exists("sNt")) list(sN=sN, long=long, sNt=sNt, outM=outM, group=group) else list(sN=sN, long=long, outM=outM, group=group)    
    return(dp) 
  }

##MA plots
maplot = function(M, A, sN, numArrays)
  {
    section = 1
    figure = 1
    sec1text = sprintf("<hr><h2><a name = \"S1\">Section %s: Individual array quality</h2></a>", section)

    app = 4 + 2*(sum(numArrays>c(4,6)))
    nfig = ceiling(numArrays/8)
    
    plotNames = paste("MA", 1:nfig, sep="")
    mapdf = paste(plotNames, "pdf", sep=".")
    mapng = paste(plotNames, "png", sep=".")
    xlimMA = quantile(A, probs=1e-4*c(1,-1)+c(0,1))
    ylimMA = quantile(M, probs=1e-4*c(1,-1)+c(0,1))
    
    dummy.df = data.frame(sN = factor(sN, levels = sN),
      x = seq_along(sN),
      y = seq_along(sN))
    
    trobj = xyplot(y ~ x | sN, dummy.df,
      xlim = xlimMA,
      ylim = ylimMA,
      xlab = "A",
      ylab = "M",
      
      panel = function(x, y, ...) {
        x <- A[, x]
        y <- M[, y]
        panel.smoothScatter(x, y, ...)
      },
      
      layout = c(app/2, 2, 1))
    id.firstpage <- seq_len(app)

    for(i in seq_len(nfig))
      {
        pdf(mapdf[i], width = 8, height = 5)
        id.thispage <- (i-1) * app + id.firstpage
        id.thispage <- id.thispage[id.thispage <= numArrays]
        print(update(trobj, index.cond = list(id.thispage)))
        dev.off()
      }
   
    png(mapng, width = 600, height = 300)
    id.thispage <- (1-1) * app + id.firstpage
    id.thispage <- id.thispage[id.thispage <= numArrays]
    print(update(trobj, index.cond = list(id.thispage)))
    dev.off()
    matext1 = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S1.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</b></center></td><td>\n", "MA plots",  basename(mapdf[1]),  basename(mapng[1]), figure)
    
    mapl = list(section=section, figure=figure, matext1=matext1, sec1text=sec1text,  nfig=nfig, mapdf=mapdf)
    return(mapl)
  }

##multidensity
multi = function(type, x, xlim, title1, title2, title3, ...)
  {
    if(type == "density")
      multidensity(x,
                   xlim = xlim,
                   main = "",
                   xlab = "",
                   ylab = "", ...)
    if(type == "ecdf")
      multiecdf(x,
                xlim = xlim,
                main = "",
                xlab = "",
                ylab = "", ...)
    mtext(title1, side = 2, adj = 0.5, padj = -4 , cex = 0.7)
    mtext(title2, side = 1, adj = 0.5, padj = 4 , cex = 0.7)
    mtext(title3, side = 3, adj = 0.5, padj = -1 , cex = 0.7)
  }

##Mapping of probes
probesmap = function(expressionset, numArrays, section, figure, dat, sN, xlim)
  {
    if(!"GC" %in% rownames(featureData(expressionset)@varMetadata))
      {
        section = section + 1
        sec3text = sprintf("<hr><h2><a name = \"S3\">Section %s: Array platform quality</a></h2>", section)          
      }
            
    figure = figure + 1

    probemapping = expressionset@featureData$hasTarget
    facgene = as.vector(probemapping)          
    facgene[probemapping == "TRUE"] = 1
    facgene[probemapping == "FALSE"] = 0
                
    gopng = "overall_GenesMapping.png"
    gopdf = "overall_GenesMapping.pdf"
          
    cols = brewer.pal(9, "Set1")
    
    png(file = gopng)
    multi("density",dat~facgene,xlim,"","","", col = cols[c(9,2)], cex.axis = 0.9)
    legend("topright", c("Mapped","Unmapped"),lty=1,lwd=2,col= c(cols[c(2,9)]), bty = "n")
    dev.copy(pdf, file = gopdf)
    dev.off()
    dev.off()

    gotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b><td><center><a name = \"S3.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</b></center></tr></td></table>\n", "Gene mapping ", basename(gopdf), basename(gopng), figure)

    promap = if(!"GC" %in% rownames(featureData(expressionset)@varMetadata)) list(section=section, figure=figure, gotext=gotext, sec3text=sec3text) else list(section=section, figure=figure, gotext=gotext) 
    return(promap)   
  }

##Heatmap
hmap = function(expressionset, sN, section, figure, outM)
  {
    section = section + 1
    sec4text = sprintf("<hr><h2><a name = \"S4\">Section %s: Between array comparison</a></h2>", section)
    figure = figure + 1
    colourRange = rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))

    d.row = as.dendrogram(hclust(outM))
    od.row = order.dendrogram(d.row)
    m = as.matrix(outM)
    rownames(m) = sN
    colnames(m) = sN

    hpdf = "heatmap.pdf"
    hpng = "heatmap.png"
    
   if("Covariate" %in% names(phenoData(expressionset)@data))
      {
        colourCov = brewer.pal(12,"Set3")
        covar = phenoData(expressionset)$Covariate
        lev = levels(as.factor(covar))
        corres = matrix(0,nrow=length(lev),ncol=2)

        png(file = hpng, w = 8*72, h = 8*72)
        print(levelplot(m[od.row,od.row],
                  scales=list(x=list(rot=90)),
                  legend=list(
                    top=list(fun=dendrogramGrob,args=list(x=d.row,side="top")),
                    right=list(fun=dendrogramGrob,args=list(x=d.row,side="right", size.add= 1, add = list(rect = list(col = "transparent", fill = colourCov[as.factor(test)])), type = "rectangle"))),
                  colorkey = list(space ="left"),
                  xlab="",ylab="",
                  col.regions=colourRange))
        
        x=0.06
        y=0.98
        
        for(i in 1:length(lev))
          {
            corres[i,] = c(unique(covar[covar == lev[i]]),colourCov[i])
            grid.text(corres[i,1], x=x, y=y, just="left")
            grid.rect(gp=gpar(fill=corres[i,2],col="transparent"), x=x-0.02, y=y, width=0.02, height=0.02)
            y=y-0.03
          }
        dev.copy(pdf, file = hpdf)
        dev.off()
        dev.off()
      } else {
        png(file = hpng, w = 8*72, h = 8*72)
        print(levelplot(m[od.row,od.row],
                  scales=list(x=list(rot=90)),
                  legend=list(
                    top=list(fun=dendrogramGrob,args=list(x=d.row,side="top")),
                    right=list(fun=dendrogramGrob,args=list(x=d.row,side="right"))),
                  colorkey = list(space ="left"),
                  xlab="",ylab="",
                  col.regions=colourRange))
        dev.copy(pdf, file = hpdf)
        dev.off()
        dev.off()
      }   
    htmltext4=sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td></table>\n", title="Heatmap representation of the distance between experiments", basename(hpdf), basename(hpng), fig = figure)           

    legendheatmap = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows a false color heatmap of between arrays distances, computed as the MAD of the M-value for each pair of arrays. <br><center><i>d<sub>ij</sub> = c &bull median|x<sub>mi</sub>-x<sub>mj</sub>|</i></center><br><br> Here, <i>x<sub>mi</sub></i> is the intensity value of the <i>m</i>-th probe on the <i>i</i>-th array, on the original data scale. <br><i>c = 1.4826</i> is a constant factor that ensures consistency with the empirical variance for Normally distributed data (see manual page of the mad function in R). This plot can serve to detect outlier arrays. <br>Consider the following decomposition of <i>x<sub>mi</sub>: x<sub>mi</sub> = z<sub>m</sub> + &beta<sub>mi</sub> + &epsilon<sub>mi</sub></i>, where <i>z<sub>m</sub></i> is the probe effect for probe <i>m</i> (the same across all arrays), <i>&epsilon<sub>mi</sub></i> are i.i.d. random variables with mean zero and <i>&beta<sub>mi</sub></i> is such that for any array <i>i</i>, the majority of values <i>&beta<sub>mi</sub></i> are negligibly small (i. e. close to zero). <i>&beta<sub>mi</sub></i> represents differential expression effects. In this model, all values <i>d<sub>ij</sub></i> are (in expectation) the same, namely <sqrt>2</sqrt> times the standard deviation of <i>&epsilon<sub>mi</i></sub> . Arrays whose distance matrix entries are way different give cause for suspicion. The dendrogram on this plot also can serve to check if, without any probe filtering, the experiments cluster accordingly to a biological meaning.</DIV></table>",  figure)
    
    hmap = list(section=section, figure=figure, htmltext4=htmltext4,sec4text=sec4text, legendheatmap=legendheatmap)
    return(hmap) 
  }

##meanSdplot
msdp = function(section, figure, con, dat)
  {
    section = section + 1
    sec5text = sprintf("<hr><h2><a name = \"S5\">Section %s: Variance mean dependency</a></h2>", section)
    figure = figure + 1
      
    mplot5 = makePlot(con=con, name = "meanSd",
      w=8, h=8, fun = function() {
        meanSdPlot(dat, cex.axis = 0.9) 
      }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td></table>\n", title="Standard deviation versus rank of the mean", fig = figure)

    htmltext5 = mplot5[[2]]
    
    legendsdmean = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\">For each feature, this plot on <b>Figure %s</b> shows the empirical standard deviation on the <i>y</i>-axis versus the rank of the mean on the <i>x</i>-axis. The red dots, connected by lines, show the running median of the standard deviation. It should be approximately horizontal, that is, show no substantial trend.</DIV>",  figure)
    
     msd = list(section=section, figure=figure, htmltext5=htmltext5,sec5text=sec5text, legendsdmean=legendsdmean)
    return(msd)   
  }

##PM.MM
pmmm = function(x, xlim, title1, title2, title3, ...)
  {
    cols = brewer.pal(9, "Set1")

    plot(density(log2(mm(x))),
         xlim = xlim,
         main = "",
         xlab = "",
         ylab = "",
         col = cols[9],...)
    lines(density(log2(pm(x))),
          xlim = xlim,
          main = "",
          xlab = "",
          ylab = "",
          col = cols[2], ...)    
    mtext(title1, side = 2, adj = 0.5, padj = -4 , cex = 0.7)
    mtext(title2, side = 1, adj = 0.5, padj = 4 , cex = 0.7)
    mtext(title3, side = 3, adj = 0.5, padj = -1 , cex = 0.7)
  }

##writing the report
 report = function(expressionset, arg, sNt, sN, sec1text, mapdf, matext1, nfig, legendMA, batext, nfig2, bapng, ftext, nfig3, fpng, legendlocal, sec2text, htmltext2, legendhom1, group, htmltext3, dtext, legendhom2, sec3text, gctext, legendgc, gotext, legendgo, sec4text, htmltext4, legendheatmap, sec5text, htmltext5, legendsdmean)
  {
### Title
    title = paste(arg$expressionset, " quality metrics report", sep="")
    titletext = sprintf("<hr><h1><center>%s</h1></center><table border = \"0\" cellspacing = 5 cellpadding = 2 >", title)
    con = openHtmlPage("QMreport", title)
    writeLines(titletext, con)

### Experiments names
    if(!is.null(sNt))
      {
        writeLines("<hr><h2>Experiment Names</h2>", con)
        for(i in seq_len(length(sN)))
          {
            if(i %in% seq(1,length(sN),by = 2))
              {
                col1 = "#e0e0f0"
                col2 = "#d0d0ff"
              }
            if(i %in% seq(2,length(sN),by = 2))
              {
                col1 = "#d0d0ff"
                col2 = "#e0e0f0"
              }
            stext = sprintf("<tr><td BGCOLOR=\"%s\"><b>%s</b></td><td BGCOLOR=\"%s\">%s</td></tr>", col1, sNt[i,2], col2, sNt[i,1])
            writeLines(stext, con)
          }
        writeLines("</table>", con)
      }

### Index
    
    writeLines("<hr><h2>Index</h2><table border = \"0\" cellspacing = 5 cellpadding = 2><UL>", con)
            
    writeLines("<tr><td><LI><b><a href=\"#S1\">Individual array quality</b></a><UL><LI><a href=\"#S1.1\">MAplot</b></a>", con)

    if("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata))
      {
        if("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))
          {
            writeLines("<LI><a href=\"#S1.2\">Spatial distribution of local background intensites</b></a>", con)
          }
        writeLines("<LI><a href=\"#S1.3\">Spatial distribution of features intensites</b></a>", con)
      }
    if(is(expressionset, "AffyBatch"))
      writeLines("<LI><a href=\"#S1.3\">Spatial distribution of features intensites</b></a>", con)

    writeLines("</UL>", con)

    writeLines( "<LI><b><a href=\"#S2\">Homogeneity between experiments</b></a><UL><LI><a href=\"#S2.1\">Boxplots</b></a><LI><a href=\"#S2.2\">Density plots</b></a></UL>", con)
            
    if("GC" %in% rownames(featureData(expressionset)@varMetadata) || "hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
      {
        writeLines("<LI><b><a href=\"#S3\">Array platform quality</b></a><UL>", con)
        if("GC" %in% rownames(featureData(expressionset)@varMetadata))
          writeLines("<LI><a href=\"#S3.1\">GC content effect</b></a>", con)
        if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
          writeLines("<LI><a href=\"#S3.2\">Gene mapping</b></a>", con)
              
        writeLines("</UL>", con)
      }
            
    writeLines("<LI><b><a href=\"#S4\">Between array comparison</b></a>" , con)
    writeLines("<LI><b><a href=\"#S5\">Variance mean dependency</b></a>" , con)

    if(is(expressionset, "AffyBatch"))
      writeLines("<LI><b><a href=\"#S6\">Affymetrix specific plots</b></a><UL><LI><a href=\"#S6.1\">RNA degradation plot</b></a><LI><a href=\"#S6.2\">RLE</b></a><LI><a href=\"#S6.3\">NUSE</b></a><LI><a href=\"#S6.4\">Affy QC stats</b></a><LI><a href=\"#S6.5\">PM MM plot</b></a></UL>" , con)

 ### Section 1
           
    writeLines("</td></tr></UL></table><table border = \"0\" cellspacing = 5 cellpadding = 2>", con)
    writeLines(sec1text, con)
    writeLines(matext1, con)
    if(nfig >= 2)
      {
        for(i in 1:nfig)
          {
            matext2 = sprintf("<A HREF=\"%s\">%s%s</A><BR>\n", basename(mapdf[i]), "MvA plot ", i)                    
            writeLines(matext2, con)
          }
      }
    writeLines("</td></table>", con)
    writeLines(legendMA, con)         
            
    if(is(expressionset, "NChannelSet") || is(expressionset, "ExpressionSet"))
      {
        if("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata))
          {
            if("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))
              {
                writeLines(batext, con)
                if(nfig2 >= 2)
                  {
                    for(i in 2:nfig2)
                      {
                        batext2 = sprintf("<A HREF=\"%s\">%s%s</A><BR>\n", basename(bapng[i]), "Spatial plots ", i)                    
                        writeLines(batext2, con)
                      }
                  }
                writeLines("</td></table>", con)
              }              
            writeLines(ftext, con)
            if(nfig3 >= 2)
              {
                for(i in 1:nfig3)
                  {
                    ftext2 = sprintf("<A HREF=\"%s\">%s%s</A><BR>\n", basename(fpng[i]), "Spatial plots ", i)
                    
                    writeLines(ftext2, con)
                  }
              }
            writeLines("</td></table>", con)
            writeLines(legendlocal, con)
          }
      }
    if(is(expressionset, "AffyBatch"))
      {
        writeLines(ftext, con)
        if(nfig3 >= 2)
          {
            for(i in 1:nfig3)
              {
                ftext2 = sprintf("<A HREF=\"%s\">%s%s</A><BR>\n", basename(fpng[i]), "Spatial plots ", i)
                    
                writeLines(ftext2, con)
              }
          }
        writeLines("</td></table>", con)
        writeLines(legendlocal, con)
      }
    
 ### Section 2
            
    writeLines(sec2text, con)
    writeLines(htmltext2, con)
    writeLines(legendhom1, con)            
    if(max(group) == 1)
      writeLines(htmltext3, con)           
    if(max(group) > 1)
      writeLines(dtext, con)              
    writeLines(legendhom2, con)
    
 ### Section 3
    if("GC" %in% rownames(featureData(expressionset)@varMetadata))
      {
        writeLines(sec3text, con)
        writeLines(gctext, con)
        writeLines(legendgc, con)
      }
    if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
      {
        if(!"GC" %in% rownames(featureData(expressionset)@varMetadata))
          writeLines(sec3text, con)          
        writeLines(gotext, con)
        writeLines(legendgo, con)
      }
            
### Section 4
    writeLines(sec4text, con)
    writeLines(htmltext4, con)
    writeLines(legendheatmap, con)
    
### Section 5
    writeLines(sec5text, con)
    writeLines(htmltext5, con)
    writeLines(legendsdmean, con)

    return(con)
  }

#################################################################
#################################################################
####################### NChannelSet Method ######################
#################################################################
#################################################################

setMethod("arrayQualityMetrics",signature(expressionset = "NChannelSet"),
          function(expressionset, outdir, force, do.logtransform, split.plots)
          {
            olddir = getwd()
            dircreation(outdir, force)

            ##data preparation
            rc = if(do.logtransform) log2(assayData(expressionset)$R) else assayData(expressionset)$R
            gc = if(do.logtransform) log2(assayData(expressionset)$G) else assayData(expressionset)$G
           
            if("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))
              {
                rbg = assayData(expressionset)$Rb
                gbg = assayData(expressionset)$Gb
              }
            ##lists to use multidensity, multiecdf and boxplot
            lredc = mat2list(rc)
            lgreenc = mat2list(gc)
      
            sN = sampleNames(expressionset)$R
            gN = featureNames(expressionset)
            numArrays = ncol(rc)
      
            dat = matrix(0, ncol = numArrays, nrow = nrow(rc))
            ##building log(ratio)
            for(a in 1:numArrays)
              dat[,a] = rc[,a] - gc[,a]
            colnames(dat) = sN
            if("dyeswap" %in% names(phenoData(expressionset)@data))
              {
                lev = levels(expressionset@phenoData$dyeswap)
                if(length(lev) != 2)
                  {
                    cat("The dyeswap slot of the phenoData must be binary.\n")
                    stop()
                  }
                reverseddye = names(expressionset@phenoData$dyeswap[expressionset@phenoData$dyeswap == min(lev)])
                dat[,reverseddye] = - dat[,reverseddye]
              }
            ldat = lapply(1:ncol(dat), function(i) dat[,i])            
     
            ##second part of data preparation
            dprep = prepdata(sN, dat, numArrays, split.plots)
            
            long = dprep$long
            sN = dprep$sN
            if("sNt" %in% names(dprep)) sNt = dprep$sNt
            outM = dprep$outM
            group = dprep$group

####################################
###Section 1 : Per array quality ###
####################################       
           
#############################
###Section 1.1 : MA-plots ###
#############################            
            ##MA-plots
            ##function from affyQCReport
            M = rc - gc
            A = 0.5*(rc + gc)
            
            MAplot = maplot(M, A, sN, numArrays)
            
            section = as.numeric(MAplot$section)
            figure = as.numeric(MAplot$figure)
            matext1 = as.character(MAplot$matext1)
            sec1text = as.character(MAplot$sec1text)
            nfig = as.numeric(MAplot$nfig)
            mapdf = as.character(MAplot$mapdf)
           
            legendMA = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> represents MA plot for each array. <br> MA plots are useful for pairwise comparisons between arrays. M and A are defined as :<br>
M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>
A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>))<br>
where I<sub>1</sub> and I<sub>2</sub> are the vectors of intensities of two channels. Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A.
Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).</DIV>", figure)          
            
#################################
###Section 1.2 : Spatial plots###
#################################

            ##Background rank representation
      
            if("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata))
              {
                r = featureData(expressionset)$X
                c = featureData(expressionset)$Y

                maxc = max(as.numeric(c))
                maxr = max(as.numeric(r))
                colourRamp = rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))
        
                if("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))
                  {
                    nfig2 = ceiling(numArrays/3)      
                    bapng = paste("background", 1:nfig2, ".png", sep="")
                    fignu = 1
                    b = 1
    
                    aR = if(maxr>maxc) maxr/maxc else maxc/maxr

                    intr = cbind(r,c,rbg)
                    intg = cbind(r,c,gbg)

                    for(a in seq_len(numArrays))
                      {
                        re = matrix(NA,ncol=maxc,nrow=maxr)
                        g = matrix(NA,ncol=maxc,nrow=maxr)
                    
                        for(i in 1:nrow(intr))
                          {
                            re[intr[i,1],intr[i,2]] = intr[i,(2+a)]
                            g[intg[i,1],intg[i,2]] = intg[i,(2+a)]
                          }
                    
                        mr = matrix(rank(re),ncol=maxc,nrow=maxr)
                        mg = matrix(rank(g),ncol=maxc,nrow=maxr)
                    
                        if(maxr>maxc){
                          mr = t(mr)
                          mg = t(mg)
                        }

                        if(a %in% seq(1,numArrays,by=3))
                          {
                            png(bapng[fignu], width = 350, height = 350*aR/2)
                            nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = FALSE),respect = FALSE)
                          }
                
                        par(xaxt = "n", yaxt = "n",mar=c(1,2,2,1))
                        image(mr, col = colourRamp)
                        mtext(sN[a],side = 3,adj = 0.5, padj = -1 ,cex = 0.7)
                        mtext("Red Intensity",side = 2, line = 0.5 ,cex = 0.7)
                        image(mg, col = colourRamp)

                        mtext("Green Intensity",side = 2, line = 0.5 ,cex = 0.7)
                        
                        if((a%%3==0) || (a == numArrays))                   
                          {
                            dev.off()
                            fignu = fignu +1
                          }
                      }
                                       
                    m = matrix(pretty(mr,9),nrow=1,ncol=length(pretty(mr,9)))
                    llbpng = "localisationlegendbackground.png"
                    png(file= llbpng, width = 3*72, height = 7*72)
                    nf <- layout(1, widths = 0.9, heights = 3, respect = TRUE)
                    image(m,xaxt="n",yaxt="n",ylab="Rank", col = colourRamp, cex.lab = 0.8, mgp = c(1.5,1,0) )
                    axis(2, label= as.list(pretty(mr,9)),at=seq(0,1,by=(1/(length(pretty(mr,9))-1))), cex.axis = 0.7, padj = 1)
                    dev.off()
               
                    batext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S1.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A></CENTER></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td><td>\n", "Spatial distribution of local background intensites", basename(bapng[1]), basename(bapng[1]), basename(llbpng))
                  }

                ##Foreground rank representation
          
                figure = figure +1
                nfig3 = ceiling(numArrays/3)      
                fpng = paste("foreground", 1:nfig3, ".png", sep="")
                fignu = 1
                b = 1
    
                aR = if(maxr>maxc) maxr/maxc else maxc/maxr

                intrf = cbind(r,c,rc)
                intgf = cbind(r,c,gc)
                for(a in seq_len(numArrays))
                  {
                    rf = matrix(NA,ncol=maxc,nrow=maxr)
                    gf = matrix(NA,ncol=maxc,nrow=maxr)
                    
                    for(i in 1:nrow(intrf))
                      {
                        rf[intrf[i,1],intrf[i,2]] = intrf[i,(2+a)]
                        gf[intgf[i,1],intgf[i,2]] = intgf[i,(2+a)]
                      }
                    
                    mrf = matrix(rank(rf),ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
                    mgf = matrix(rank(gf),ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))

                    if(maxr>maxc){
                      mrf = t(mrf)
                      mgf = t(mgf)
                    }

                    if(a %in% seq(1,numArrays,by=3))
                      {
                        png(fpng[fignu], width = 350, height = 350*aR/2)
                        nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = FALSE),
                                     respect = FALSE)
                      }

                    par(xaxt = "n", yaxt = "n",mar=c(1,2,2,1))
                    image(mrf, col = colourRamp)
                    mtext(sN[a],side = 3,adj = 0.5, padj = -1 ,cex = 0.7)
                    mtext("Red Intensity",side = 2, line = 0.5 ,cex = 0.7)
                    image(mgf, col = colourRamp)

                    mtext("Green Intensity",side = 2, line = 0.5 ,cex = 0.7)
                    if((a%%3==0) || (a == numArrays))                   
                      {
                        dev.off()
                        fignu = fignu +1
                      }
                  }
                              
                m = matrix(pretty(mrf,9),nrow=1,ncol=length(pretty(mrf,9)))
                llfpng = "localisationlegendforeground.png"
                png(file= llfpng, width = 3*72, height = 7*72)
                nf <- layout(1, widths = 0.9, heights = 3, respect = TRUE)

                image(m,xaxt="n",yaxt="n",ylab="Rank", col = colourRamp, cex.lab = 0.8, mgp = c(1.5,1,0) )
                axis(2, label= as.list(pretty(mrf,9)),at=seq(0,1,by=(1/(length(pretty(mrf,9))-1))), cex.axis = 0.7, padj = 1)
                dev.off()
                
                ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name = \"S1.3\"><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></a></A><center><b>Figure %s</b></center></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td><td>\n", "Spatial distribution of features intensites", basename(fpng[1]), basename(fpng[1]), figure, basename(llfpng))

                legendlocal = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s:</b> False color representations of the arrays' spatial distributions of feature intensities and, if available, local background estimates. The color scale is shown in the panel on the right, and it is proportional to the ranks. These plots may help in identifying patterns that may be caused, for example, spatial gradients in the hybridization chamber, air bubbles, spotting or plating problems.</DIV>", figure)          
              }

#################################################
###Section 2 : Homogeneity between experiments###
#################################################

#############################
###Section 2.1 : Boxplots ###
#############################
            
            section = section + 1
            sec2text = sprintf("<hr><h2><a name = \"S2\">Section %s: Homogeneity between experiments</a></h2>", section)
            
            figure = figure + 1
            
            xlim = c(min(na.omit(dat)),max(na.omit(dat)))
            xlimr = c(min(na.omit(rc)),max(na.omit(rc)))
            xlimg = c(min(na.omit(gc)),max(na.omit(gc)))

            colours = brewer.pal(12, "Paired")
            if(numArrays <= 50)
              {
                xname = sN
                mai = c(0,0.4,0.2,0.2)
                xaxt = "s"
              }
            if(numArrays > 50)
              {
                xname = FALSE
                mai = c(0,0.4,0.2,0.2)
                xaxt = "n"
              }

            mplot2 = makePlot(con=con, name = "boxplot",
                     w=15, h=8, fun = function() {
                       nf = layout(matrix(1:3,1,3,byrow=TRUE),
                                    c(2,2,2), 2, TRUE)
                       par(cex.axis = 1,
                           pty = "s",
                           lheight = ((1/log10(numArrays))*long),
                           mai = mai,
                           omi = c(0,0,0,0),
                           xaxt = xaxt)
                       boxplot(lgreenc, col = colours[4], las = 3, range = 0,
                               names = xname, ylim = c(min(c(rc,gc)),max(c(rc,gc))), title = "Green Channel")
                       boxplot(lredc, col = colours[6], las = 3, range = 0,
                               names = xname, ylim =  c(min(c(rc,gc)),max(c(rc,gc))), title = "Red Channel")
                       boxplot(ldat, col = colours[2], las = 3, range = 0,
                               names = xname, ylim = xlim, title = "Log(Ratio)")
                     }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</center></b></td></table>\n", title="Boxplots", fig = figure)

            htmltext2 = mplot2[[2]]
            legendhom1 = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> presents boxplots. On the left panel, the green boxes correspond to the log<sub>2</sub> intensities of the green channel. On the middle panel the red boxes correspond to the log<sub>2</sub> intensities of the red channel. The right panel shows the boxplots of log<sub>2</sub>(ratio).</DIV>", figure)          

      
############################
###Section 2.2 : Density ###
############################
            
            ##Density if 1 group
            figure = figure + 1

            if(max(group) == 1)
              {
                mplot3 = makePlot(con=con, name = "density",
                         w=10, h=10, fun = function() {
                           nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),
                                        c(2.3,2,2), c(2,2), TRUE)
                           par(mar=c(0,5,1,1),xaxt = "n", cex.axis = 0.9)
                           multi("density",lgreenc,xlimg,"Density","","Green Channel")
                           par(mar=c(0,2,1,1),xaxt = "n")
                           multi("density",lredc,xlimr,"","","Red Channel")
                           par(mar=c(0,2,1,1),xaxt = "n")
                           multi("density",ldat,xlim,"","","Log(Ratio)")
                           par(mar=c(1,5,0,1), xaxt = "s")
                           multi("ecdf",lgreenc,xlimg,"ECDF","log(intensity)","")
                           par(mar=c(1,2,0,1), xaxt = "s")
                           multi("ecdf",lredc,xlimr,"","log(intensity)","")
                           par(mar=c(1,2,0,1), xaxt = "s")
                           multi("ecdf",ldat,xlim,"","log(ratio)","")}, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", title="Density plots", fig = figure)
                htmltext3 = mplot3[[2]]          
              }
      
            ##Density if more than 1 group
            if(max(group) > 1)
              {
                dpng = "density.png"
                dpdf = "density.pdf"
                pdf(file = dpdf)
                for(n in 1:max(group))
                  {
                    nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),
                                 c(2.3,2,2), c(2,2), TRUE)
              
                    par(mar=c(0,5,1,1),xaxt = "n", cex.axis = 0.9)
                    multi("density",lgreenc[group==n],xlimg,"Density","","Green Channel")
                    par(mar=c(0,2,1,1),xaxt = "n")
                    multi("density",lredc[group==n],xlimr,"","","Red Channel")
                    mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -3 ,cex = 1)
                    par(mar=c(0,2,1,1),xaxt = "n")
                    multi("density",ldat[group==n],xlim,"","","Log(Ratio)")
                    par(mar=c(1,5,0,1), xaxt = "s")
                    multi("ecdf",lgreenc[group==n],xlimg,"ECDF","log(intensity)","")
                    par(mar=c(1,2,0,1), xaxt = "s")
                    multi("ecdf",lredc[group==n],xlimr,"","log(intensity)","")
                    par(mar=c(1,2,0,1), xaxt = "s")
                    multi("ecdf",ldat[group==n],xlim,"","log(ratio)","")
                  }
                dev.off()
                png(file = dpng)
                nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),
                             c(2.3,2,2), c(2,2), TRUE)
              
                par(mar=c(0,5,1,1),xaxt = "n", cex.axis = 0.9)
                multi("density",lgreenc[group==1],xlimg,"Density","","Green Channel")
                par(mar=c(0,2,1,1),xaxt = "n")
                multi("density",lredc[group==1],xlimr,"","","Red Channel")
                mtext("Group 1",side = 3,adj = 0.5, padj = -3 ,cex = 1)
                par(mar=c(0,2,1,1),xaxt = "n")
                multi("density",ldat[group==1],xlim,"","","Log(Ratio)")
                par(mar=c(1,5,0,1), xaxt = "s")
                multi("ecdf",lgreenc[group==1],xlimg,"ECDF","log(intensity)","")
                par(mar=c(1,2,0,1), xaxt = "s")
                multi("ecdf",lredc[group==1],xlimr,"","log(intensity)","")
                par(mar=c(1,2,0,1), xaxt = "s")
                multi("ecdf",ldat[group==1],xlim,"","log(ratio)","")
                dev.off()
                dtext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", "Density plots", basename(dpdf), basename(dpng), figure)
              }
            legendhom2 = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows density estimates (histograms) of the data. Arrays whose distributions are very different from the others should be considered for possible problems.</DIV>", figure)

########################################
###Section 3 : Array platform quality###
########################################
     
#######################
###Section 3.1 : GC ###
#######################
    
            if("GC" %in% rownames(featureData(expressionset)@varMetadata))
              {
                section = section + 1
                sec3text = sprintf("<hr><h2><a name = \"S3\">Section %s: Array platform quality</a></h2>", section)    
            
                figure = figure + 1
                gcpdf = "GCcontent.pdf"
                pdf(file = gcpdf)
                ngc = c(2:9)
                colb = brewer.pal(9, "Blues")
                colg = brewer.pal(9, "Greens")
                colr = brewer.pal(9, "OrRd")         
                fac = round(as.numeric(as.matrix(featureData(expressionset)$GC)),-1)
         
                nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9),3,3,byrow = FALSE),
                             c(1.2,1.2,1.2), c(2,1.8,2), FALSE)
          
                for ( a in 1:numArrays )
                  {
                    par(mar=c(0,3.5,3.5,1))
                    multi("density",rc[,a]~fac,xlimr,"Density","","Red Channel", col = colr[ngc],  xaxt = "n")
                    par(mar=c(2,3.5,0,1))
                    multi("ecdf",rc[,a]~fac,xlimr,"ECDF","","", col = colr[ngc])
                    par(mar=c(2,3.5,2,1))
                    boxplot(rc[,a]~fac, col = colr[ngc], range =  0, main ="")
                    mtext("Boxplot",side = 2,adj = 0.5, padj = -4 ,cex = 0.7)
              
                    par(mar=c(0,3.5,3.5,1))
                    multi("density",gc[,a]~fac,xlimg,"","","Green Channel", col = colg[ngc],  xaxt = "n")
                    mtext(sN[a],side = 3,adj = 0.5, padj = -2 ,cex = 1, font = 2)
                    par(mar=c(2,3.5,0,1))
                    multi("ecdf",gc[,a]~fac,xlimg,"","","", col = colg[ngc])
                    par(mar=c(2,3.5,2,1))
                    boxplot(gc[,a]~fac, col = colg[ngc], lwd = 1, range =  0, main = "")
                    par(mar=c(0,3.5,3.5,1))
                    multi("density",dat[,a]~fac,xlim,"","","Log(ratio)", col = colb[ngc],  xaxt = "n")
                    par(mar=c(2,3.5,0,1))
                    multi("ecdf",dat[,a]~fac,xlim,"","","", col = colb[ngc])
                    par(mar=c(2,3.5,2,1))
                    boxplot(dat[,a]~fac, col = colb[ngc], lwd = 1, range = 0, main = "")
                  }
                dev.off()
        
                gcopng = "overall_GCcontent.png"
                gcopdf = "overall_GCcontent.pdf"
                png(file = gcopng)
                nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9),3,3,byrow = FALSE),
                             c(1.2,1.2,1.2), c(2,1.8,2), FALSE)
                par(mar=c(0,3.5,3.5,1))
                multi("density",rc~fac,xlimr,"Density","","Red Channel", col = colr[ngc],  xaxt = "n")
                par(mar=c(2,3.5,0,1))
                multi("ecdf",rc~fac,xlimr,"ECDF","","", col = colr[ngc])
                par(mar=c(2,3.5,2,1))
                boxplot(rc~fac, col = colr[ngc], range =  0, main ="")
                mtext("Boxplot",side = 2,adj = 0.5, padj = -4 ,cex = 0.7)
            
                par(mar=c(0,3.5,3.5,1))
                multi("density",gc~fac,xlimg,"","","Green Channel", col = colg[ngc],  xaxt = "n")
                par(mar=c(2,3.5,0,1))
                multi("ecdf",gc~fac,xlimg,"","","", col = colg[ngc])
                par(mar=c(2,3.5,2,1))
                boxplot(gc~fac, col = colg[ngc], lwd = 1, range =  0, main = "")
            
                par(mar=c(0,3.5,3.5,1))
                multi("density",dat~fac,xlim,"","","Log(ratio)", col = colb[ngc],  xaxt = "n")
                par(mar=c(2,3.5,0,1))
                multi("ecdf",dat~fac,xlim,"","","", col = colb[ngc])
                par(mar=c(2,3.5,2,1))
                boxplot(dat~fac, col = colb[ngc], lwd = 1, range = 0, main = "")
                dev.copy(pdf, file = gcopdf)
                dev.off()
                dev.off()
          
                gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><center><a name = \"S3.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</b></center></tr></td></table>\n", "GC content effect ", basename(gcpdf), "per array", " and ", basename(gcopdf), "global", basename(gcopdf), basename(gcopng), figure)
                legendgc = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the distributions of the log<sub>2</sub> intensities grouped by the percentage of cytosines (C) and guanines (G) among the nucleotides in each probe. From the top to the bottom, kernel density estimates, empirical cumulative distribution functions (ECDF) and boxplots are represented. Box and line colors in the three panels correspond to the same groups. Cytosine and guanine are able to form three hydrogen bonds, while adenine (A) and thymine (T) only form two, hence oligonucleotides with a higher proportion of C and G can form more stable hybridization bindings. This should result in higher intensities measured on the array, regardless of the abundance of target molecules.</DIV>",  figure)
              }

######################################
###Section 3.1 : Mapping of probes ###
######################################
    
            if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
              {
                pmap = probesmap(expressionset, numArrays, section, figure, dat, sN, xlim)

                section = pmap$section
                figure = pmap$figure
                gotext = pmap$gotext
                if(!"GC" %in% rownames(featureData(expressionset)@varMetadata)) sec3text = as.character(pmap$sec3text)                
                legendgo = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the density distributions of the log<sub>2</sub> ratios grouped by the mapping of the probes. Blue, density estimate of log<sub>2</sub> ratios of probes annotated \"TRUE\" in the <b>\"hasTarget\"</b> slot. Gray, probes annotated \"FALSE\" in the <b>\"hasTarget\"</b> slot.</DIV>",  figure)    
              }
            
##########################################
###Section 4 : Between array comparison###
##########################################
            
            ##Heatmap
            hmap = hmap(expressionset, sN, section, figure, outM)
            
            section = hmap$section
            figure = hmap$figure
            htmltext4 = hmap$htmltext4
            sec4text = hmap$sec4text
            legendheatmap = hmap$legendheatmap
            
##########################################
###Section 5 : Variance Mean Dependency###
##########################################
            
            ##meanSdplot
            msdp = msdp(section, figure, con, dat)
            
            section = msdp$section
            figure = msdp$figure
            htmltext5 = msdp$htmltext5
            sec5text = msdp$sec5text
            legendsdmean = msdp$legendsdmean

##########################
### Writing the report ###
##########################

            arg = as.list(match.call(expand.dots = TRUE))
            con = if(exists("sNt")) report(expressionset, arg, sNt, sN, sec1text, mapdf, matext1, nfig, legendMA, batext, nfig2, bapng, ftext, nfig3, fpng, legendlocal, sec2text, htmltext2, legendhom1, group, htmltext3, dtext, legendhom2, sec3text, gctext, legendgc, gotext, legendgo, sec4text, htmltext4, legendheatmap, sec5text, htmltext5, legendsdmean) else report(expressionset, arg, sN, sNt = NULL, sec1text, mapdf, matext1, nfig, legendMA, batext, nfig2, bapng, ftext, nfig3, fpng, legendlocal, sec2text, htmltext2, legendhom1, group, htmltext3, dtext, legendhom2, sec3text, gctext, legendgc, gotext, legendgo, sec4text, htmltext4, legendheatmap, sec5text, htmltext5, legendsdmean)

            writeLines("</table>", con)
            closeHtmlPage(con)
            setwd(olddir)
            
          }####end set method NChannelSet
          )

#################################################################
#################################################################
#################### ExpressionSet Functions ####################
#################################################################
#################################################################

aqm.expressionset = function(expressionset, outdir = getwd(), force = FALSE, do.logtransform = FALSE, split.plots = FALSE, arg)
  {
    olddir = getwd()
    dircreation(outdir, force)
   
    ##data preparation
    dat = if(do.logtransform) log2(exprs(expressionset)) else exprs(expressionset)
   
    sN = sampleNames(expressionset)
    ##list to use multidensity, multiecdf and boxplot
    ldat = mat2list(dat)
    gN = featureNames(expressionset)
    numArrays = ncol(dat)
    
    ##second part of data preparation
    dprep = prepdata(sN, dat, numArrays, split.plots)
            
    long = dprep$long
    sN = dprep$sN
    if("sNt" %in% names(dprep)) sNt = dprep$sNt
    outM = dprep$outM
    group = dprep$group
      
####################################
###Section 1 : Per array quality ###
####################################
    
#############################
###Section 1.1 : MA-plots ###
#############################
    
    ##MA-plots
    ##function from affyQCReport   
    medArray = rowMedians(na.omit(dat))
    M =  na.omit(dat)-medArray
    A =  (na.omit(dat)+medArray)/2

    MAplot = maplot(M, A, sN, numArrays)

    section = as.numeric(MAplot$section)
    figure = as.numeric(MAplot$figure)
    matext1 = as.character(MAplot$matext1)
    sec1text = as.character(MAplot$sec1text)
    nfig = as.numeric(MAplot$nfig)
    mapdf = as.character(MAplot$mapdf)   
    
    legendMA = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> represents MA plot for each array. <br> MA plots are useful for pairwise comparisons between arrays. M and A are defined as :<br>
M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>
A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>))<br>
where I<sub>1</sub> and I<sub>2</sub> are the vectors of intensities of two channels. Rather than comparing each array to every other array, here we compare each array to a single median  \"pseudo\"-array. Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A. Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).</DIV>", figure)     

#################################
###Section 1.2 : Spatial Plots###
#################################
  
    if(is(expressionset, "AffyBatch"))
      {
        figure = figure + 1
        maxc = ncol(expressionset)
        maxr = nrow(expressionset)
        nfig3 = ceiling(numArrays/6)      
        colourRamp <- rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))

        fpng = paste("foreground", 1:nfig3, ".png", sep="")
        fignu = 1
        aR = if(maxr>maxc) maxr/maxc else maxc/maxr

        for(a in seq_len(numArrays))
          {
            rfi = rank(dat[,a])
            mrfi = matrix(rfi,ncol=maxc,nrow=maxr,byrow=T)                                        
            if(maxr>maxc)
              mrfi = t(mrfi)
                
            if(a %in% seq(1,numArrays,by=6))
              {
                png(fpng[fignu], width = 500, height = 500*aR/2)
                nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = TRUE),
                             widths = c(1.5,1.5,1.5),
                             heights = c(1.5*aR,1.5*aR),
                             respect = TRUE)
              }
                
            par(xaxt = "n", yaxt = "n",mar=c(1,2,2,1))
            image(mrfi, col = colourRamp)
            mtext(sN[a],side = 3, adj = 0.5, padj = -1 ,cex = 0.7)
                
            if((a%%6==0) || (a == numArrays))                   
              {
                dev.off()
                fignu = fignu +1
              }
          }
                                   
        m = matrix(pretty(mrfi,9),nrow=1,ncol=length(pretty(mrfi,9)))
        llpng = "localisationlegend.png" 
        png(file= llpng, width = 3*72, height = 7*72)
        nf <- layout(1, widths = 0.9, heights = 3, respect = TRUE)
        image(m,xaxt="n",yaxt="n",ylab="Rank", col = colourRamp, cex.lab = 0.8, mgp = c(1.5,1,0) )
        axis(2, label= as.list(pretty(mrfi,9)),at=seq(0,1,by=(1/(length(pretty(mrfi,9))-1))), cex.axis = 0.7, padj = 1)
        dev.off()
           
        ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name=\"S1.2\"><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></a></A><center><b>Figure %s</b></center></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td><td>\n", "Spatial distribution of features intensities", basename(fpng[1]), basename(fpng[1]), figure, basename(llpng))

        legendlocal = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s:</b> False color representations of the arrays' spatial distributions of feature intensities. The color scale is shown in the panel on the right, and it is proportional to the ranks. These plots may help in identifying patterns that may be caused, for example, spatial gradients in the hybridization chamber, air bubbles, spotting or plating problems.</DIV>", figure)
      }

    if(is(expressionset, "ExpressionSet"))
      {
        if("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata))
          {
            r = featureData(expressionset)$X
            c = featureData(expressionset)$Y

            maxc = max(as.numeric(c))
            maxr = max(as.numeric(r))
            colourRamp = rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))
          
            figure = figure +1
            nfig3 = ceiling(numArrays/6)      
            fpng = paste("foreground", 1:nfig3, ".png", sep="")
            fignu = 1
            b = 1
    
            aR = if(maxr>maxc) maxr/maxc else maxc/maxr

            intf = cbind(r,c,dat)
            for(a in seq_len(numArrays))
              {
                fg = matrix(NA,ncol=maxc,nrow=maxr)
                   
                for(i in 1:nrow(intf))
                  {
                    fg[intf[i,1],intf[i,2]] = intf[i,(2+a)]
                  }
                    
                mfg = matrix(rank(fg),ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
                    
                if(maxr>maxc){
                  mfg = t(mfg)
                }

                if(a %in% seq(1,numArrays,by=6))
                  {
                    png(fpng[fignu], width = 350, height = 350*aR/2)
                    nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = FALSE),
                                 respect = FALSE)
                  }

                par(xaxt = "n", yaxt = "n",mar=c(1,2,2,1))
                image(mfg, col = colourRamp)
                mtext(sN[a],side = 3, line = 0.5 ,cex = 0.7)
                if((a%%6==0) || (a == numArrays))                   
                  {
                    dev.off()
                    fignu = fignu +1
                  }
              }
                              
            m = matrix(pretty(mfg,9),nrow=1,ncol=length(pretty(mfg,9)))
            llfpng = "localisationlegendforeground.png"
            png(file= llfpng, width = 3*72, height = 7*72)
            nf <- layout(1, widths = 0.9, heights = 3, respect = TRUE)

            image(m,xaxt="n",yaxt="n",ylab="Rank", col = colourRamp, cex.lab = 0.8, mgp = c(1.5,1,0) )
            axis(2, label= as.list(pretty(mfg,9)),at=seq(0,1,by=(1/(length(pretty(mfg,9))-1))), cex.axis = 0.7, padj = 1)
            dev.off()
                
            ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name = \"S1.3\"><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></a></A><center><b>Figure %s</b></center></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td><td>\n", "Spatial distribution of features intensites", basename(fpng[1]), basename(fpng[1]), figure, basename(llfpng))

            legendlocal = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s:</b> False color representations of the arrays' spatial distributions of feature intensities. The color scale is shown in the panel on the right, and it is proportional to the ranks. These plots may help in identifying patterns that may be caused, for example, spatial gradients in the hybridization chamber, air bubbles, spotting or plating problems.</DIV>", figure)

          }
      }

#################################################
###Section 2 : Homogeneity between experiments###
#################################################

#############################
###Section 2.1 : Boxplots ###
#############################

    section = section + 1
    sec2text = sprintf("<hr><h2><a name = \"S2\">Section %s: Homogeneity between experiments</h2></a>", section)
            
    figure = figure + 1
    mplot2 = makePlot(con=con, name = "boxplot",
      w=8, h=8, fun = function() {
        colours = brewer.pal(12, "Paired")
        par(cex.axis = 1, pty = "s", lheight =((1/log10(numArrays))*long), mai = c(((long/12)+0.2),0.4,0.2,0.2) , omi = c(0,0,0,0))
        boxplot(ldat, col = colours[2], las = 3, range = 0, names = sN)
      }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</center></b></td></table>\n", title = "Boxplots", fig = figure)
    htmltext2 = mplot2[[2]]
    legendhom1 = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> presents boxplots of the data.</DIV>", figure)          
    
############################
###Section 2.2 : Density ###
############################
    
    ##Density if 1 group
    figure = figure + 1
    xlim = c(min(na.omit(dat)),max(na.omit(dat)))
    if(max(group) == 1)
      {
        mplot3 = makePlot(con=con, name = "density",
          w=10, h=10, fun = function() {
            xlim = c(min(na.omit(dat)),max(na.omit(dat)))
            nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(2.8,2.8),c(1.8,2), TRUE)
            par(xaxt = "n", cex.axis = 0.8, mar = c(0,5,2,5))
            multi("density",ldat,xlim,"Density","","")
            par(xaxt = "s", cex.axis = 0.8, mar = c(4,5,0,5))
            multi("ecdf",ldat,xlim,"ECDF","log(intensity)","")
          }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", title="Density plots", fig = figure)
        htmltext3 = mplot3[[2]]
      }
                
    ##Density if more than 1 group
    if(max(group) > 1)
      {
        dpng = "density.png"
        dpdf = "density.pdf"
        pdf(file = dpdf)
        xlim = c(min(na.omit(dat)),max(na.omit(dat)))
        for(n in 1:max(group))
          {
            nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(2.8,2.8),c(1.8,2), TRUE)
            par(xaxt = "n", cex.axis = 0.8, mar = c(0,5,2,5))
            multi("density",ldat[group==n],xlim,"Density","","")
            mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -1 ,cex = 1)
            par(xaxt = "s", cex.axis = 0.8, mar = c(4,5,0,5))
            multi("ecdf",ldat[group==n],xlim,"ECDF","log(intensity)","")
          }
        dev.off()
        png(file = dpng)
        nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(2.8,2.8),c(1.8,2), TRUE)
        par(xaxt = "n", cex.axis = 0.8, mar = c(0,5,2,5))
        multi("density",ldat[group==1],xlim,"Density","","")
        mtext("Group 1",side = 3,adj = 0.5, padj = -1 ,cex = 1)
        par(xaxt = "s", cex.axis = 0.8, mar = c(4,5,0,5))
        multi("ecdf",ldat[group==1],xlim,"ECDF","log(intensity)","")
        dev.off()
        dtext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></CENTER><BR></tr></td></table>\n", "Density plots", basename(dpdf), basename(dpng), figure)
      }
    legendhom2 = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows density estimates (histograms) of the data. Arrays whose distributions are very different from the others should be considered for possible problems.</DIV>", figure)          

########################################
###Section 3 : Array platform quality###
########################################

#######################
###Section 3.1 : GC ###
#######################
    
    if("GC" %in% rownames(featureData(expressionset)@varMetadata))
      {
        section = section + 1
        sec3text = sprintf("<hr><h2><a name = \"S3\">Section %s: Array platform quality</h2></a>", section)    
        figure = figure + 1
        gcpdf = "GCcontent.pdf"
        pdf(file = gcpdf)
        ngc = c(2:9)
        colb = brewer.pal(9, "Blues")
        colg = brewer.pal(9, "Greens")
        colr = brewer.pal(9, "OrRd")         
        fac = round(as.numeric(as.matrix(featureData(expressionset)$GC)),-1)
        
        nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9),3,3,byrow = FALSE),
                     c(1.2,1.2,1.2), c(2,1.8,2), FALSE)
        
        for ( a in 1:numArrays )
          {
            par(mar=c(0,3.5,3.5,1))
            multi("density",dat[,a]~fac,xlim,"Density","",sN[a], col = colb[ngc],  xaxt = "n")
            par(mar=c(2,3.5,0,1))
            multi("ecdf",dat[,a]~fac,xlim,"ECDF","","", col = colb[ngc],  xaxt = "n")
            par(mar=c(2,3.5,2,1))
            boxplot(dat[,a]~fac, col = colb[ngc], lwd = 1, range = 0)
            mtext("Boxplot",side = 2,adj = 0.5, padj = -3.3 ,cex = 0.7)
          }
        dev.off()
        
        gcopng = "overall_GCcontent.png"
        gcopdf = "overall_GCcontent.pdf"
        png(file = gcopng)
        nf <- layout(matrix(c(1,2,3),3,1,byrow = FALSE),
                     2.2, c(2,1.8,2), TRUE)
        par(mar=c(0,3.5,3.5,1))
        multi("density",dat~fac,xlim,"Density","","Log(ratio)", col = colb[ngc],  xaxt = "n")
        par(mar=c(2,3.5,0,1))
        multi("ecdf",dat~fac,xlim,"ECDF","","", col = colb[ngc])
        par(mar=c(2,3.5,2,1))
        boxplot(dat~fac, col = colb[ngc], lwd = 1, range = 0)
        mtext("Boxplot",side = 2,adj = 0.5, padj = -3.3 ,cex = 0.7)
        dev.copy(pdf, file = gcopdf)
        dev.off()
        dev.off()
                
        gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><center><a name = \"S3.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</b></center></tr></td></table>\n", "GC content effect ", basename(gcpdf), "per array", " and ", basename(gcopdf), "global", basename(gcopdf), basename(gcopng), figure)
        legendgc = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the distributions of the log<sub>2</sub> intensities grouped by the percentage of cytosines (C) and guanines (G) among the nucleotides in each probe. From the top to the bottom, kernel density estimates, empirical cumulative distribution functions (ECDF) and boxplots are represented. Box and line colors in the three panels correspond to the same groups. Cytosine and guanine are able to form three hydrogen bonds, while adenine (A) and thymine (T) only form two, hence oligonucleotides with a higher proportion of C and G can form more stable hybridization bindings. This should result in higher intensities measured on the array, regardless of the abundance of target molecules.</DIV>",  figure)
      }

######################################
###Section 3.1 : Mapping of probes ###
######################################
    
    if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
      {
        pmap = probesmap(expressionset, numArrays, section, figure, dat, sN, xlim)

        section = pmap$section
        figure = pmap$figure
        gotext = pmap$gotext
        sec3text = as.character(pmap$sec3text)
                
        legendgo = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the density distributions of the log<sub>2</sub> intensities grouped by the mapping of the probes. Blue, density estimate of intensities of probes annotated \"TRUE\" in the <b>\"hasTarget\"</b> slot. Gray, probes annotated \"FALSE\" in the <b>\"hasTarget\"</b> slot.</DIV>", figure)          
      }
    
##########################################
###Section 4 : Between array comparison###
##########################################
   
    ##Heatmap
    hmap = hmap(expressionset, sN, section, figure, outM)            
    section = hmap$section
    figure = hmap$figure
    htmltext4 = hmap$htmltext4
    sec4text = hmap$sec4text
    legendheatmap = hmap$legendheatmap

##########################################
###Section 5 : Variance Mean Dependency###
##########################################
    
    ##meanSdplot
    msdp = msdp(section, figure, con, dat)            
    section = msdp$section
    figure = msdp$figure
    htmltext5 = msdp$htmltext5
    sec5text = msdp$sec5text
    legendsdmean = msdp$legendsdmean    

#########################
### Writing the report###
#########################
    
    con = if(exists("sNt")) report(expressionset, arg, sNt, sN, sec1text, mapdf, matext1, nfig, legendMA, batext, nfig2, bapng, ftext, nfig3, fpng, legendlocal, sec2text, htmltext2, legendhom1, group, htmltext3, dtext, legendhom2, sec3text, gctext, legendgc, gotext, legendgo, sec4text, htmltext4, legendheatmap, sec5text, htmltext5, legendsdmean) else report(expressionset, arg, sN, sNt = NULL, sec1text, mapdf, matext1, nfig, legendMA, batext, nfig2, bapng, ftext, nfig3, fpng, legendlocal, sec2text, htmltext2, legendhom1, group, htmltext3, dtext, legendhom2, sec3text, gctext, legendgc, gotext, legendgo, sec4text, htmltext4, legendheatmap, sec5text, htmltext5, legendsdmean)


    l = list(numArrays=numArrays,sN=sN, section=section, figure=figure, con=con, dat=dat, olddir=olddir)

    return(l) 
  }

#################################################################
#################################################################
###################### ExpressionSet Method #####################
#################################################################
#################################################################

setMethod("arrayQualityMetrics",signature(expressionset="ExpressionSet"),
          function(expressionset, outdir, force, do.logtransform, split.plots)
          {
            arg = as.list(match.call(expand.dots = TRUE))

            l = aqm.expressionset(expressionset, outdir, force, do.logtransform, split.plots, arg)
            con = l[[5]]
            olddir = l[[7]]
            
            writeLines("</table>", con)
            closeHtmlPage(con)
            setwd(olddir)
          })##end set method ExpressionSet

#################################################################
#################################################################
######################## AffyBatch Method #######################
#################################################################
#################################################################

setMethod("arrayQualityMetrics",signature(expressionset="AffyBatch"),
          function(expressionset, outdir, force, do.logtransform,
          split.plots){
            
            arg = as.list(match.call(expand.dots = TRUE))

            l = aqm.expressionset(expressionset, outdir, force, do.logtransform, split.plots, arg)
            numArrays = as.numeric(l$numArrays)
            sN = l$sN
            section = as.numeric(l$section)
            figure = as.numeric(l$figure)
            con = l$con
            dat = l$dat
            olddir = l$olddir              

############################
###Section 8 : Affy plots###
############################

            cols = brewer.pal(9, "Set1")            
            section = section + 1
            sec6text = sprintf("<hr><h2><a name = \"S6\">Section %s: Affymetrix specific plots</h2></a>", section)
            writeLines(sec6text, con)
            
            figure1 = figure + 1
            acol = sample(brewer.pal(8, "Dark2"), numArrays, replace = (8<numArrays))
            rnaDeg = AffyRNAdeg(expressionset)
            affypng1 = "RNAdeg.png"
            affypdf1 = "RNAdeg.pdf"
            png(file = affypng1)
            plotAffyRNAdeg(rnaDeg, cols = acol, lwd = 2)
            dev.copy(pdf, file = affypdf1)
            dev.off()
            dev.off()
                
            figure2 = figure1 + 1
            pp1 = preprocess(expressionset)
            dataPLM = fitPLM(pp1, background = FALSE, normalize = FALSE)
            affypng2 = "RLE.png"
            affypdf2 = "RLE.pdf"
            png(file = affypng2)
            Mbox(dataPLM, ylim = c(-1, 1), names = sN, col = cols[2],
                 whisklty = 0, staplelty = 0, main = "RLE", las = 3, cex.axis = 0.8)
            dev.copy(pdf, file = affypdf2)
            dev.off()
            dev.off()
            
            figure3 = figure2 + 1
            affypng3 = "NUSE.png"
            affypdf3 = "NUSE.pdf"
            png(file = affypng3)
            boxplot(dataPLM, ylim = c(0.95, 1.5), names = sN,
                    outline = FALSE, col = cols[2], main = "NUSE", las = 2, cex.axis = 0.8)
            dev.copy(pdf, file = affypdf3)
            dev.off()
            dev.off()
            
            figure4 = figure3 + 1
            affypng4 = "qc.png"
            affypdf4 = "qc.pdf"
            png(file = affypng4)
            plot(qc(expressionset), cex.axis = 0.8)
            dev.copy(pdf, file = affypdf4)
            dev.off()
            dev.off()

            affytext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name = \"S6.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><center><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><a name = \"S6.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><center><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><a name = \"S6.3\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><center><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><a name = \"S6.4\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><center><BR><b>Figure %s</b></CENTER></tr></td></table>\n", "RNA degradation plot", basename(affypdf1), basename(affypng1), figure1, "RLE plot", basename(affypdf2), basename(affypng2), figure2, "NUSE plot", basename(affypdf3), basename(affypng3), figure3, "Diagnostic plot recommended by Affy", basename(affypdf4), basename(affypng4), figure4)
            writeLines(affytext, con)
            legendaffy = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\">In this section we present diagnostic plots based on tools provided in the affyPLM package. In <b>Figure %s</b> a RNA digestion plot is computed on normalized data (so that standard deviation is equal to 1). In this plot each array is represented by a single line. It is important to identify any array(s) that has a slope which is very different from the others. The indication is that the RNA used for that array has potentially been handled quite differently from the other arrays.  <b>Figure %s</b> is a Relative Log Expression (RLE) plot and an array that has problems will either have larger spread, or will not be centered at M = 0, or both. <b>Figure %s</b> is a Normalized Unscaled Standard Error (NUSE) plot. Low quality arrays are those that are significantly elevated or more spread out, relative to the other arrays. NUSE values are not comparable across data sets. Both RLE and NUSE are performed on preprocessed data (background correction and quantile normalization). <b>Figure %s</b> represents the diagnostic plot recommended by Affymetrix. It is fully describe in the simpleaffy.pdf vignette of the package simpleaffy. Any metrics (circles and triangles) that is shown in red is out of the manufacturer's specific boundaries and suggests a potential problem.</DIV>", figure1, figure2, figure3, figure4)            
            writeLines(legendaffy, con)
            
            ##PM.MM
            figure5 = figure4 + 1            
            cols = brewer.pal(9, "Set1")
            xlim = c(min(na.omit(dat)),max(na.omit(dat)))            
            pmopng = "overall_PM.MM.png"
            pmopdf = "overall_PM.MM.pdf"       
        
            png(file = pmopng)
            pmmm(expressionset,xlim,"","","", cex.axis = 0.9)
            legend("topright", c("PM","MM"),lty=1,lwd=2,col= c(cols[c(2,9)]), bty = "n")
            dev.copy(pdf, file = pmopdf)
            dev.off()
            dev.off()
            
            pmotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b><td><center><a name = \"S6.5\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><br><b>Figure %s</b></center></tr></td></table>\n", "Perfect matchs and mismatchs ", basename(pmopdf), basename(pmopng), figure5)
            legendpmo = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the density distributions of the log<sub>2</sub> intensities grouped by the matching of the probes. Blue, density estimate of intensities of perfect match probes and gray the mismatch probes.</DIV>",  figure5)
            writeLines(pmotext, con)
            writeLines(legendpmo, con)

            writeLines("</table>", con)
            closeHtmlPage(con)
            setwd(olddir)

          })##end set method AffyBatch
