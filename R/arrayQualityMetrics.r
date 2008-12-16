setGeneric("arrayQualityMetrics",
           function(expressionset,
                    outdir = getwd(),
                    force = FALSE,
                    do.logtransform = FALSE,
                    split.plots = FALSE,
                    intgroup = "Covariate",
                    grouprep = FALSE)
           standardGeneric("arrayQualityMetrics"))

#################################################################
#################################################################
########################### FUNCTIONS ###########################
#################################################################
#################################################################

##lists
mat2list = function(x)
      lapply(seq_len(ncol(x)), function(i) x[,i])

##log
logtransform = function(x)
  ifelse(x>0, suppressWarnings(log2(x)), rep(NA_real_, length(x)))
##  {
##    xl = log2(x)
##    xl[xl == "-Inf"] = NA
##    return(xl)
##  }

##makePlot
makePlot = function(con, name, w, h, fun, title, text, fig) {
  outf = paste(name, c("pdf", "png"), sep=".")
  nrppi = 72

  pdf(outf[1], width=w, height=h, pointsize=14)
  fun()
  dev.off()

  wd = w*nrppi
  hg = h*nrppi
  
  png(outf[2], width=wd, height=hg, pointsize=14)
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
              
      if(!force && length(outdirContents)>0)
        stop(sprintf("'%s' is not empty.", outdir))
      if(force || length(outdirContents)==0)
        message(sprintf("The report will be written in '%s'. ", outdir))
        setwd(outdir)
       } else {
        dir.create(outdir, recursive=TRUE)
        message(sprintf("The directory '%s' has been created.", outdir))
        setwd(outdir)
      }
  }

##second part of data preparation
prepdata = function(sN, dat, numArrays, split.plots)
  {
    ##maximum length of the arrays names to adjust margins and font size on the plot axis
    long = max(nchar(sN))
      
    ##if arrays names are to long they are replaced by figures
    sNt = cbind(sN,seq_len(length(sN)))        
    colnames(sNt)=c("SampleName","New denomination")
    if(long >= 15 || (numArrays >= 50 && long >= 10))
      {
        sN = seq_len(length(sN))
        long = max(nchar(sN))
      }

    outM = as.dist(dist2(dat))
      
    k = if(split.plots) split.plots else k = numArrays
    ##attribute randomly the arrays to different groups
    group = sample(rep(seq_len(ceiling(numArrays/k)),k),numArrays)
    
    dp = list(sN=sN, long=long, sNt=sNt, outM=outM, group=group)
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
    
    plotNames = paste("MA", seq_len(nfig), sep="")
    mapdf = paste(plotNames, "pdf", sep=".")
    mapng = paste(plotNames, "png", sep=".")
    xlimMA = quantile(A, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
    ylimMA = quantile(M, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
    
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
        panel.smoothScatter(x, y ,...)
      },
      as.table=TRUE,      
      layout = c(app/2, 2, 1),
      strip = function(..., bg) strip.default(..., bg ="#cce6ff"))
    
    id.firstpage = seq_len(app)

    for(i in seq_len(nfig))
      {
        pdf(mapdf[i], width = 8, height = 5)
        id.thispage = (i-1) * app + id.firstpage
        id.thispage = id.thispage[id.thispage <= numArrays]
        print(update(trobj, index.cond = list(id.thispage)))
        dev.off()
      }

    png(mapng, width = 600, height = 300)
    id.thispage = id.firstpage[id.firstpage <= numArrays]
    print(update(trobj, index.cond = list(id.thispage)))
    dev.off()
    matext1 = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S1.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</b></center></td><td>\n", "MA plots",  basename(mapdf[1]),  basename(mapng[1]), figure)
    
    mapl = list(section=section, figure=figure, matext1=matext1, sec1text=sec1text,  nfig=nfig, mapdf=mapdf)
    return(mapl)
  }

##multidensity
multi = function(type, x, xlim, title1, title2, title3, ...)
  {
    switch(type,
           "density" ={
             multidensity(x,
                          xlim = xlim,
                          main = "",
                          xlab = "",
                          ylab = "", ...)},
           "ecdf"={
             multiecdf(x,
                       xlim = xlim,
                       main = "",
                       xlab = "",
                       ylab = "", ...)
             mtext(title1, side = 2, adj = 0.5, padj = -4 , cex = 0.9)
             mtext(title2, side = 1, adj = 0.5, padj = 4 , cex = 0.9)
             mtext(title3, side = 3, adj = 0.5, padj = -1 , cex = 0.9)
           })
  }

##Mapping of probes
probesmap = function(expressionset, numArrays, section, figure, dat, sN, xlim, type,con)
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

    mplot6 = makePlot(con=con, name = "ProbesMapping",
      w=8, h=8, fun = function() {
        meanSdPlot(dat, cex.axis = 0.9, ylab = "Standard deviation of the intensities", xlab="Rank(mean of intensities)")        
        cols = brewer.pal(9, "Set1")
        if(type == 2)
          xla = "log(ratio)"    
        if(type == 1)
          xla = "log(intensity)"
        multi("density",dat~facgene,xlim,"Density",xla,"", col = cols[c(9,2)], cex.axis = 0.7)
        legend("topright", legend=c("Mapped","Unmapped"),lty=1,lwd=2,col= c(cols[c(2,9)]), bty = "n", cex =0.9)}, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b><td><center><a name = \"S3.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</b></center></tr></td></table>\n",title="Probes mapping", fig = figure)

    gotext = mplot6[[2]]

    promap = if(!"GC" %in% rownames(featureData(expressionset)@varMetadata)) list(section=section, figure=figure, gotext=gotext, sec3text=sec3text) else list(section=section, figure=figure, gotext=gotext) 
    return(promap)   
  }

##Heatmap
hmap = function(expressionset, sN, section, figure, outM, numArrays, intgroup)
  {
    section = section + 1
    sec4text = sprintf("<hr><h2><a name = \"S4\">Section %s: Between array comparison</a></h2>", section)
    figure = figure + 1
    colourRange = rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))

    d.row = as.dendrogram(hclust(outM, method = "single"))
    od.row = order.dendrogram(d.row)
    m = as.matrix(outM)
    colnames(m) = sN
    rownames(m) = sN
   
    hpdf = "heatmap.pdf"
    hpng = "heatmap.png"
    
   if(all(intgroup %in% names(phenoData(expressionset)@data)))
      {
        covar = lapply(1:length(intgroup), function(i) pData(expressionset)[colnames(pData(expressionset))==intgroup[i]][,1])
        lev = lapply(1:length(intgroup), function(i) levels(as.factor(unlist(covar[i]))))
        corres = lapply(1:length(intgroup), function(i) matrix(0,nrow=length(lev[[i]]),ncol=2))
        colourCov = lapply(1:length(intgroup), function(i) brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[i]))

        key = lapply(1:length(intgroup), function(i) list(rect = list(col=unlist(colourCov[i])[as.factor(levels(as.factor(unlist(covar[i]))))]), text = list(levels(as.factor(unlist(covar[i]))))))
        key = unlist(key, recursive=F)
        key$rep = FALSE
      
        foo = draw.key(key = key)


        
        hfig = levelplot(m[od.row,od.row],
          scales=list(x=list(rot=90)),
          legend=list(
            top=list(fun=dendrogramGrob,args=list(x=d.row,side="top")),
            right=list(fun=dendrogramGrob,
              args=list(x=d.row,side="right", size.add=1, 
                add = sapply(1:length(intgroup), function(i) list(rect = list(col="transparent",fill = unlist(colourCov[i])[as.factor(unlist(covar[i]))]))),
                type = "rectangle"))),
          colorkey = list(space ="left"),
          xlab="",ylab="",
          col.regions=colourRange,
          main = foo)
        
      } else {
        hfig = levelplot(m[od.row,od.row],
          scales=list(x=list(rot=90)),
          legend=list(
            top=list(fun=dendrogramGrob,args=list(x=d.row,side="top")),
            right=list(fun=dendrogramGrob,args=list(x=d.row,side="right"))),
          colorkey = list(space ="left"),
          xlab="",ylab="",
          col.regions=colourRange)
      }
    png(file = hpng, w = 8*72, h = 8*72)
    print(hfig)
    dev.off()
    pdf(file = hpdf)
    print(hfig)
    dev.off()
    
    htmltext4=sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td></table>\n", title="Heatmap representation of the distance between arrays", basename(hpdf), basename(hpng), fig = figure)
    
    leghmspe = if(is(expressionset, "BeadLevelList")) "the values used are the summarized ones obtained by using the function createBeadSummaryData from the package beadarray." else "without preprocessing." 

    legendheatmap = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows a false color heatmap of between arrays distances, computed as the mean absolute difference (L<sub>1</sub>-distance) of the vector of M-values for each pair of arrays. <br><center><i>d<sub>xy</sub> = mean|M<sub>xi</sub>-M<sub>yi</sub>|</i></center><br><br> Here, <i>M<sub>xi</sub></i> is the M-value of the <i>i</i>-th probe on the <i>x</i>-th array, %s <br>This plot can serve to detect outlier arrays. <br>Consider the following decomposition of <i>M<sub>xi</sub>: M<sub>xi</sub> = z<sub>i</sub> + &beta;<sub>xi</sub> + &epsilon;<sub>xi</sub></i>, where <i>z<sub>i</sub></i> is the probe effect for probe <i>i</i> (the same across all arrays), <i>&epsilon;<sub>xi</sub></i> are i.i.d. random variables with mean zero and <i>&beta;<sub>xi</sub></i> is such that for any array <i>x</i>, the majority of values <i>&beta;<sub>xi</sub></i> are negligibly small (i. e. close to zero). <i>&beta;<sub>xi</sub></i> represents differential expression effects. In this model, all values <i>d<sub>xy</sub></i> are (in expectation) the same, namely <sqrt>2</sqrt> times the standard deviation of <i>&epsilon;<sub>xi</i></sub> . Arrays whose distance matrix entries are way different give cause for suspicion. The dendrogram on this plot also can serve to check if, without any probe filtering, the arrays cluster accordingly to a biological meaning. The colour scale is chosen to cover the range of L1-distances encountered in the dataset.</DIV></table>",  figure, leghmspe)
    
    if(all(intgroup %in% names(phenoData(expressionset)@data)))
      {
        hmap = list(section=section, figure=figure, htmltext4=htmltext4,sec4text=sec4text, legendheatmap=legendheatmap, numArrays=numArrays, intgroup=intgroup)
      } else {
        hmap = list(section=section, figure=figure, htmltext4=htmltext4,sec4text=sec4text, legendheatmap=legendheatmap, numArrays=numArrays)
      }
    return(hmap) 
  }

##meanSdplot
msdp = function(expressionset, section, figure, con, dat)
  {
    section = section + 1
    sec5text = sprintf("<hr><h2><a name = \"S5\">Section %s: Variance mean dependence</a></h2>", section)
    figure = figure + 1
      
    mplot5 = makePlot(con=con, name = "meanSd",
      w=8, h=8, fun = function() {
        meanSdPlot(dat, cex.axis = 0.9, ylab = "Standard deviation of the intensities", xlab="Rank(mean of intensities)") 
      }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td></table>\n", title="Standard deviation versus rank of the mean", fig = figure)

    htmltext5 = mplot5[[2]]

    legsdspe = if(is(expressionset, "BeadLevelList")) "For each bead type obtained by createBeadSummaryData from the package beadarray," else "For each feature,"
    
    legendsdmean = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\">%s the plot on <b>Figure %s</b> shows the empirical standard deviation of the intensities of all the arrays on the <i>y</i>-axis versus the rank of the mean of intensities of the arrays on the <i>x</i>-axis. The red dots, connected by lines, show the running median of the standard deviation. After vsn normalization, this should be approximately horizontal, that is, show no substantial trend.</DIV>", legsdspe, figure)
    
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
         xlab = "log(Intensity)",
         ylab = "Density",
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

##boxplot.stats2 is the same as boxplot.stats except that the boxes that are less wide than the others are not detected as outliers
##  FIXME - unclear what "boxes ... are not detected as outliers" means. 
boxplot.stats2 = function (x, coef = 1.5, do.conf = TRUE, do.out = TRUE) 
{
  if (coef < 0) 
    stop("'coef' must not be negative")
  nna <- !is.na(x)
  n <- sum(nna)
  stats <- stats::fivenum(x, na.rm = TRUE)
  iqr <- diff(stats[c(2, 4)])
  if (coef == 0) 
    do.out <- FALSE
  else {
    out <- if (!is.na(iqr)) {
      ##x < (stats[2] - coef * iqr) |
      x > (stats[4] + coef *iqr)
    }
    else !is.finite(x)
    if (any(out[nna], na.rm = TRUE)) 
      stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
  }
  conf <- if (do.conf) 
    stats[3] + c(-1.58, 1.58) * iqr/sqrt(n)
  list(stats = stats, n = n, conf = conf, out = if (do.out) x[out & nna] else numeric(0))
}

##SCORES COMPUTATION
scores = function(expressionset, numArrays, M, ldat, outM, dat, rc, gc, maxc, maxr, nuse, rle)
  {
    ## MA plot
    mamean = colMeans(abs(M), na.rm=TRUE)
    mastat = boxplot.stats(mamean)
    maout = sapply(seq_len(length(mastat$out)), function(x) which(mamean == mastat$out[x]))
    ## boxplot
    b = if(is(expressionset,"BeadLevelList")) {
      boxplotBeads(expressionset, plot = FALSE,
                   whatToPlot = c(single="G", two="M")[expressionset@arrayInfo$channels])
    } else {
      boxplot(ldat, plot = FALSE, range = 0)
    }
    bmeanstat = boxplot.stats(b$stat[3,])
    bmeanout = sapply(seq_len(length(bmeanstat$out)), function(x) which(b$stat[3,] == bmeanstat$out[x]))

    biqr = b$stat[4,] -  b$stat[2,] 
    biqrstat = boxplot.stats(biqr)
    biqrout = sapply(seq_len(length(biqrstat$out)), function(x) which(biqr == biqrstat$out[x]))

    ## heatmap
    madsum  = rowSums(as.matrix(outM), na.rm=TRUE)
    madstat = boxplot.stats2(madsum)
    madout = sapply(seq_len(length(madstat$out)), function(x) which(madsum == madstat$out[x]))

    if(!is(expressionset, "BeadLevelList") && (is(expressionset, "AffyBatch") || ("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata))))
      {
        ## spatial plot
        if(is(expressionset, "AffyBatch") || is(expressionset, "ExpressionSet"))
          {
            mdat = lapply(seq_len(numArrays), function(x) matrix(dat[,x],ncol=maxc,nrow=maxr,byrow=T))
            loc = sapply(mdat, function(x) {
              apg = abs(fft(x)) ## absolute values of the periodogramme
              lowFreq = apg[1:4, 1:4]
              lowFreq[1,1] = 0  # drop the constant component
              highFreq = c(apg[-(1:4), ], apg[1:4, -(1:4)])
              return(sum(lowFreq)/sum(highFreq))
            })
            locstat = boxplot.stats(loc)
            locout = sapply(seq_len(length(locstat$out)), function(x) which(loc == locstat$out[x]))
          }
        if(is(expressionset, "NChannelSet"))
          {
            mdatr = lapply(seq_len(numArrays), function(x) matrix(rc[,x],ncol=maxc,nrow=maxr,byrow=T))
            locr = sapply(mdatr, function(x) {
              apg = abs(fft(x)) ## absolute values of the periodogramme
              lowFreq = apg[1:4, 1:4]
              lowFreq[1,1] = 0  # drop the constant component
              highFreq = c(apg[-(1:4), ], apg[1:4, -(1:4)])
              return(sum(lowFreq)/sum(highFreq))
            })
            locstatr = boxplot.stats(locr)
            locoutr = sapply(seq_len(length(locstatr$out)), function(x) which(locr == locstatr$out[x]))
            
            mdatg = lapply(seq_len(numArrays), function(x) matrix(gc[,x],ncol=maxc,nrow=maxr,byrow=T))
            locg = sapply(mdatg, function(x) {
              apg = abs(fft(x)) ## absolute values of the periodogramme
              lowFreq = apg[1:4, 1:4]
              lowFreq[1,1] = 0  # drop the constant component
              highFreq = c(apg[-(1:4), ], apg[1:4, -(1:4)])
              return(sum(lowFreq)/sum(highFreq))
            })
            locstatg = boxplot.stats(locg)
            locoutg = sapply(seq_len(length(locstatg$out)), function(x) which(locg == locstatg$out[x]))
          }
      }

    if(is(expressionset, "AffyBatch"))
      {
        ## RLE
        rlemed = rle$stats[3,]
        rleout = which(abs(rlemed) > 0.1)

        ## NUSE
        nusemeanstat = boxplot.stats(nuse$stat[3,])
        nusemeanout = sapply(seq_len(length(nusemeanstat$out)), function(x) which(nuse$stat[3,] == nusemeanstat$out[x]))

        nuseiqr = nuse$stat[4,] -  nuse$stat[2,] 
        nuseiqrstat = boxplot.stats(nuseiqr)
        nuseiqrout = sapply(seq_len(length(nuseiqrstat$out)), function(x) which(nuseiqr == nuseiqrstat$out[x]))


        scoresout = list(maout,locout,union(bmeanout,biqrout),madout,rleout,union(nusemeanout,nuseiqrout))
      }
    
    if(!is(expressionset, "BeadLevelList") && ("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata)) && is(expressionset, "ExpressionSet") )      
      scoresout = list(maout,locout,union(bmeanout,biqrout),madout)
    
    if(!is(expressionset, "BeadLevelList") && ("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata)) && is(expressionset, "NChannelSet") )      
      scoresout = list(maout,union(locoutr,locoutg),union(bmeanout,biqrout),madout)
    
    if(is(expressionset, "BeadLevelList") || (!is(expressionset, "AffyBatch") && (!("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata)))))
      scoresout = list(maout,union(bmeanout,biqrout),madout)
          
    return(scoresout)
  }

##writing the report
report = function(expressionset, arg, sNt, sN, sec1text, mapdf, matext1, nfig, legendMA, batext, nfig2, bapng, ftext, pttext, legendpt, nfig3, fpng, legendlocal, sec2text, htmltext2, legendhom1, group, htmltext3, legendhom2, sec3text, gctext, legendgc, gotext, legendgo, sec4text, htmltext4, legendheatmap, sec5text, htmltext5, legendsdmean, scores)
  {
### Title
        argum = arg$expressionset
        title = paste(deparse(substitute(argum)), " quality metrics report", sep="")
    titletext = sprintf("<hr><h1><center>%s</h1></center><table border = \"0\" cellspacing = 5 cellpadding = 2 style=\"font-size: 13; font-family: Lucida Grande; text-align:center\">", title)
    con = openHtmlPage("QMreport", title)
    writeLines(titletext, con)

### Arrays names
    writeLines("<hr><h2>Summary</h2>", con)
    col1 = "#d0d0ff"
    col2 = "#e0e0f0"
    if(is(expressionset, "AffyBatch"))
      stext1 = sprintf("<tr BGCOLOR=\"%s\"><td><b>Array &#35;</b></td><td><b>Array Name</b></td><td><b>MA-plot</b></td><td><b>Spatial distribution</b></td><td><b>Boxplots/Density plots</b></td><td><b>Heatmap</b></td><td><b>RLE</b></td><td><b>NUSE</b></td></tr>", col1)

    if(!is(expressionset, "BeadLevelList") && ("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata)))
      stext1 = sprintf("<tr BGCOLOR=\"%s\"><td><b>Array &#35;</b></td><td><b>Array Name</b></td><td><b>MA-plot</b></td><td><b>Spatial distribution</b></td><td><b>Boxplots/Density plots</b></td><td><b>Heatmap</b></td></tr>", col1)


    if(is(expressionset, "BeadLevelList") || (!is(expressionset, "AffyBatch") && (!("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata)))))
      stext1 = sprintf("<tr BGCOLOR=\"%s\"><td><b>Array &#35;</b></td><td><b>Array Name</b></td><td><b>MA-plot</b></td><td><b>Boxplots/Density plots</b></td><td><b>Heatmap</b></td></tr>", col1)
    
    writeLines(stext1, con)
    for(i in seq_len(length(sN)))
      {
        if(i %in% seq(1,length(sN),by = 2))
          col = col2
        if(i %in% seq(2,length(sN),by = 2))
          col = col1
        if(is(expressionset, "AffyBatch"))
          stext2 = sprintf("<tr BGCOLOR=\"%s\"><td><b>%s</b></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>", col, sNt[i,2], sNt[i,1], scores[i,1], scores[i,2], scores[i,3], scores[i,4], scores[i,5], scores[i,6])
        
        if(!is(expressionset, "BeadLevelList") && ("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata)))
          stext2 = sprintf("<tr BGCOLOR=\"%s\"><td><b>%s</b></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>", col, sNt[i,2],sNt[i,1], scores[i,1], scores[i,2], scores[i,3], scores[i,4])       
        
        if(is(expressionset, "BeadLevelList") || (!is(expressionset, "AffyBatch") && (!("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata)))))
          stext2 = sprintf("<tr BGCOLOR=\"%s\"><td><b>%s</b></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>", col, sNt[i,2], sNt[i,1], scores[i,1], scores[i,2], scores[i,3])
        
        writeLines(stext2, con)
      }
    writeLines("</table>", con)
    
    scoreslegend =  sprintf("*array identified as having a potential problem or as being an outlier.")
    
    writeLines(scoreslegend, con)
      
    
### Index
    
    writeLines("<hr><h2>Index</h2><table border = \"0\" cellspacing = 5 cellpadding = 2><UL>", con)
            
    writeLines("<tr><td><LI><b><a href=\"#S1\">Individual array quality</b></a><UL><LI><a href=\"#S1.1\">MAplot</b></a>", con)

    if(!is(expressionset, "BeadLevelList") && "X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata))
      {
        if("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))
          writeLines("<LI><a href=\"#S1.2\">Spatial distribution of local background intensites</b></a>", con)
        
        writeLines("<LI><a href=\"#S1.3\">Spatial distribution of feature intensites</b></a>", con)
        writeLines("<LI><a href=\"#S1.4\">Row-Column effect</b></a>", con)
      }
    
    if(is(expressionset, "AffyBatch"))
      {
        writeLines("<LI><a href=\"#S1.2\">Spatial distribution of feature intensites</b></a>", con)
        maxc = ncol(expressionset)
        maxr = nrow(expressionset)       
        if(maxr*maxc < 1000000)
          writeLines("<LI><a href=\"#S1.4\">Row-Column effect</b></a>", con)
      }

    if(is(expressionset, "BeadLevelList"))
      {
        writeLines("<LI><a href=\"#S1.2\">Spatial distribution of local background intensites</b></a>", con)
        writeLines("<LI><a href=\"#S1.3\">Spatial distribution of feature intensites</b></a>", con)
      }

    writeLines("</UL>", con)

    writeLines( "<LI><b><a href=\"#S2\">Homogeneity between arrays</b></a><UL><LI><a href=\"#S2.1\">Boxplots</b></a><LI><a href=\"#S2.2\">Density plots</b></a></UL>", con)
            
    if((!is(expressionset, "BeadLevelList")) && ("GC" %in% rownames(featureData(expressionset)@varMetadata) || "hasTarget" %in% rownames(featureData(expressionset)@varMetadata)))
      {
        writeLines("<LI><b><a href=\"#S3\">Array platform quality</b></a><UL>", con)
        if("GC" %in% rownames(featureData(expressionset)@varMetadata))
          writeLines("<LI><a href=\"#S3.1\">GC content effect</b></a>", con)
        if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
          writeLines("<LI><a href=\"#S3.2\">Gene mapping</b></a>", con)
              
        writeLines("</UL>", con)
      }
            
    writeLines("<LI><b><a href=\"#S4\">Between array comparison</b></a>" , con)
    writeLines("<LI><b><a href=\"#S5\">Variance mean dependence</b></a>" , con)

    if(is(expressionset, "AffyBatch"))
      writeLines("<LI><b><a href=\"#S6\">Affymetrix specific plots</b></a><UL><LI><a href=\"#S6.1\">RNA degradation plot</b></a><LI><a href=\"#S6.2\">RLE</b></a><LI><a href=\"#S6.3\">NUSE</b></a><LI><a href=\"#S6.4\">Affymetrix QC stats</b></a><LI><a href=\"#S6.5\">PM MM plot</b></a></UL>" , con)

 ### Section 1
           
    writeLines("</td></tr></UL></table><table border = \"0\" cellspacing = 5 cellpadding = 2>", con)
    writeLines(sec1text, con)
    writeLines(matext1, con)
    if(nfig >= 2)
      {
        for(i in seq_len(nfig))
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
                for(i in seq_len(nfig3))
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
            for(i in seq_len(nfig3))
              {
                ftext2 = sprintf("<A HREF=\"%s\">%s%s</A><BR>\n", basename(fpng[i]), "Spatial plots ", i)
                    
                writeLines(ftext2, con)
              }
          }
        writeLines("</td></table>", con)
        writeLines(legendlocal, con)
        maxc = ncol(expressionset)
        maxr = nrow(expressionset)
        if(maxc*maxr < 1000000)
          {
            writeLines(pttext, con)
            writeLines(legendpt, con)
          }
      }
    
    if(is(expressionset, "BeadLevelList"))
      {
        writeLines(batext, con)
        if(nfig2 >= 2)
          {
            for(i in seq_len(nfig2))
              {
                batext2 = sprintf("<A HREF=\"%s\">%s%s</A><BR>\n", basename(bapng[i]), "Spatial plots ", i)                    
                writeLines(batext2, con)
              }
          }
        writeLines("</td></table>", con)
        
        writeLines(ftext, con)
        if(nfig3 >= 2)
          {
            for(i in seq_len(nfig3))
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
    writeLines(htmltext3, con)           
    writeLines(legendhom2, con)
    
 ### Section 3
    if((!is(expressionset, "BeadLevelList")) && ("GC" %in% rownames(featureData(expressionset)@varMetadata)))
      {
        writeLines(sec3text, con)
        writeLines(gctext, con)
        writeLines(legendgc, con)
      }
    if((!is(expressionset, "BeadLevelList")) && ("hasTarget" %in% rownames(featureData(expressionset)@varMetadata)))
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
########################## RGList Method ########################
#################################################################
#################################################################

setMethod("arrayQualityMetrics",signature(expressionset = "RGList"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            expressionset = try(as(expressionset, "NChannelSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a RGList and cannot be converted automatically in a NChannelSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          })####end set method RGList

#################################################################
#################################################################
######################## marrayRaw Method #######################
#################################################################
#################################################################
setMethod("arrayQualityMetrics",signature(expressionset = "marrayRaw"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            expressionset = try(as(expressionset, "NChannelSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a marrayRaw and cannot be converted automatically in a NChannelSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          })####end set method marrayRaw

#################################################################
#################################################################
######################### MAList Method #########################
#################################################################
#################################################################
setMethod("arrayQualityMetrics",signature(expressionset = "MAList"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            expressionset = try(as(expressionset, "ExpressionSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a MAList and cannot be converted automatically in a ExpressionSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          })####end set method MAList

#################################################################
#################################################################
######################## marrayNorm Method ######################
#################################################################
#################################################################
setMethod("arrayQualityMetrics",signature(expressionset = "marrayNorm"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            expressionset = try(as(expressionset, "ExpressionSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a marrayNorm and cannot be converted automatically in a ExpressionSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          })####end set method marrayNorm


#################################################################
#################################################################
####################### NChannelSet Method ######################
#################################################################
#################################################################

setMethod("arrayQualityMetrics",signature(expressionset = "NChannelSet"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            
            olddir = getwd()
            on.exit(setwd(olddir))

            ##data preparation
            if(do.logtransform)
              {
                rc = logtransform(assayData(expressionset)$R)
                gc = logtransform(assayData(expressionset)$G)
              } else {
                rc = assayData(expressionset)$R
                gc = assayData(expressionset)$G
              }
            
            dircreation(outdir, force)
            
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
      
            ##building log(ratio)
            dat = rc - gc
            colnames(dat) = sN
            if("dyeswap" %in% names(phenoData(expressionset)@data))
              {
                lev = levels(expressionset@phenoData$dyeswap)
                if(length(lev) != 2)
                  stop("The dyeswap slot of the phenoData must be binary.\n")
                reverseddye = names(expressionset@phenoData$dyeswap[expressionset@phenoData$dyeswap == min(lev)])
                dat[,reverseddye] = - dat[,reverseddye]
              }
            ldat = mat2list(dat)
     
            ##second part of data preparation
            dprep = prepdata(sN, dat, numArrays, split.plots)
            
            long = dprep$long
            sN = dprep$sN
            sNt = dprep$sNt
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
            A = 0.5*(rc+gc)
            
            MAplot = maplot(M, A, sN, numArrays)
            
            section = as.numeric(MAplot$section)
            figure = as.numeric(MAplot$figure)
            matext1 = as.character(MAplot$matext1)
            sec1text = as.character(MAplot$sec1text)
            nfig = as.numeric(MAplot$nfig)
            mapdf = as.character(MAplot$mapdf)
           
            legendMA = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the MA plot for each array. M and A are defined as:<br>M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>)),<br>where I<sub>1</sub> and I<sub>2</sub> are the vectors of intensities of the two channels. Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A. Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).</DIV>", figure)          
            
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
                    bapng = paste("background", seq_len(nfig2), ".png", sep="")
                    fignu = 1
                    b = 1
    
                    aR = if(maxr>maxc) maxr/maxc else maxc/maxr

                    intr = cbind(r,c,rbg)
                    intg = cbind(r,c,gbg)

                    for(a in seq_len(numArrays))
                      {
                        re = matrix(NA,ncol=maxc,nrow=maxr)
                        g = matrix(NA,ncol=maxc,nrow=maxr)
                    
                        for(i in seq_len(nrow(intr)))
                          {
                            re[intr[i,1],intr[i,2]] = intr[i,(2+a)]
                            g[intg[i,1],intg[i,2]] = intg[i,(2+a)]
                          }
                    
                        mr = matrix(rank(re),ncol=maxc,nrow=maxr)
                        mg = matrix(rank(g),ncol=maxc,nrow=maxr)
                    
                        if(maxr > maxc){
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
                        mtext("Rank(red intensity)",side = 2, line = 0.5 ,cex = 0.7)
                        image(mg, col = colourRamp)

                        mtext("Rank(green intensity)",side = 2, line = 0.5 ,cex = 0.7)
                        
                        if((a%%3 == 0) || (a == numArrays))                   
                          {
                            dev.off()
                            fignu = fignu +1
                          }
                      }
                                       
                    m = matrix(pretty(mr,9),nrow=1,ncol=length(pretty(mr,9)))
                    llbpng = "localisationlegendbackground.png"
                    png(file= llbpng, width = 3*72, height = 7*72)
                    nf <- layout(1, widths = 0.9, heights = 3, respect = TRUE)
                    image(m,xaxt="n",yaxt="n",ylab="Rank(Intensity)", col = colourRamp, cex.lab = 0.8)
                    axis(2, label= as.list(pretty(mr,9)),at=seq(0,1,by=(1/(length(pretty(mr,9))-1))), cex.axis = 0.7, padj = 1)
                    dev.off()
               
                    batext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S1.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A></CENTER></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td><td>\n", "Spatial distribution of local background intensites", basename(bapng[1]), basename(bapng[1]), basename(llbpng))
                  }

###################################
###Section 1.3 : Spatial plots 2###
###################################

                ##Foreground rank representation
          
                figure = figure +1
                nfig3 = ceiling(numArrays/3)      
                fpng = paste("foreground", seq_len(nfig3), ".png", sep="")
                fignu = 1
                b = 1
    
                aR = if(maxr > maxc) maxr/maxc else maxc/maxr

                intrf = cbind(r,c,rc)
                intgf = cbind(r,c,gc)
                for(a in seq_len(numArrays))
                  {
                    rf = matrix(NA,ncol=maxc,nrow=maxr)
                    gf = matrix(NA,ncol=maxc,nrow=maxr)
                    
                    for(i in seq_len(nrow(intrf)))
                      {
                        rf[intrf[i,1],intrf[i,2]] = intrf[i,(2+a)]
                        gf[intgf[i,1],intgf[i,2]] = intgf[i,(2+a)]
                      }
                    
                    mrf = matrix(rank(rf),ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
                    mgf = matrix(rank(gf),ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))

                    if(maxr > maxc){
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
                    mtext("Rank(red intensity)",side = 2, line = 0.5 ,cex = 0.7)
                    image(mgf, col = colourRamp)

                    mtext("Rank(green intensity)",side = 2, line = 0.5 ,cex = 0.7)
                    if((a%%3 == 0) || (a == numArrays))                   
                      {
                        dev.off()
                        fignu = fignu +1
                      }
                  }
                              
                m = matrix(pretty(mrf,9),nrow=1,ncol=length(pretty(mrf,9)))
                llfpng = "localisationlegendforeground.png"
                png(file= llfpng, width = 3*72, height = 7*72)
                nf <- layout(1, widths = 0.9, heights = 3, respect = TRUE)

                image(m,xaxt="n",yaxt="n",ylab="Rank(Intensity)", col = colourRamp, cex.lab = 0.8)
                axis(2, label= as.list(pretty(mrf,9)),at=seq(0,1,by=(1/(length(pretty(mrf,9))-1))), cex.axis = 0.7, padj = 1)
                dev.off()
                
                ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name = \"S1.3\"><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></a></A><center><b>Figure %s</b></center></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td><td>\n", "Spatial distribution of feature intensites", basename(fpng[1]), basename(fpng[1]), figure, basename(llfpng))

                legendlocal = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s:</b> False color representations of the arrays' spatial distributions of feature intensities and, if available, local background estimates. The color scale is shown in the panel on the right, and it is proportional to the ranks of the probe intensities. This has indeed the potential to detect patterns that are small in amplitude, but systematic within array. These may not be consequential for the downstream data analysis, but if properly interpreted, could e.g. still be useful for technology and experimental protocol optimisation as it helps in identifying patterns that may be caused by, for example, spatial gradients in the hybridization chamber, air bubbles, spotting or plating problems.</DIV>", figure)          
              
###############################
###Section 1.4 : R/C effect ###
###############################
    
                figure = figure + 1
                ptpdf = "row_column.pdf"
                ptpng = "row_column.png"
                colours = brewer.pal(12, "Paired")

                facx = round(as.numeric(as.matrix(featureData(expressionset)$X)),-1)
                facy = round(as.numeric(as.matrix(featureData(expressionset)$Y)),-1)         
                pdf(file = ptpdf, width=6, height=6)
                nf <- layout(matrix(1:6,3,2,byrow = FALSE),
                             c(2.1,2), c(2,2,2), FALSE)
          
                for ( a in seq_len(numArrays))
                  {
                    par(mar=c(2,3.5,3.5,1), cex.axis=0.7)
                    boxplot(rc[,a]~facx, col = colours[6], range =  0, main ="", whisklty = 0, staplelty = 0)
                    mtext(sN[a], side = 3,padj=-1,adj=1.1,cex=1)
                    mtext("Row", side = 3,padj=-0.5,cex=0.8)
                    mtext("log(red intensity)", side = 2, adj = 0.5,padj=-4,cex=0.7)

                    boxplot(gc[,a]~facx, col = colours[4], lwd = 1, range =  0, main = "", whisklty = 0, staplelty = 0)
                    mtext("log(green intensity", side = 2, adj = 0.5,padj=-4,cex=0.7)
                    boxplot(dat[,a]~facx, col = colours[2], lwd = 1, range = 0, main = "", whisklty = 0, staplelty = 0)
                    mtext("log(ratio)", side = 2, adj = 0.5,padj=-4,cex=0.7)

                    boxplot(rc[,a]~facy, col = colours[6], range =  0, main ="", whisklty = 0, staplelty = 0)
                    mtext("Column", side = 3,padj=-0.5,cex=0.8)

                    boxplot(gc[,a]~facy, col = colours[4], lwd = 1, range =  0, main = "", whisklty = 0, staplelty = 0)
                    boxplot(dat[,a]~facy, col = colours[2], lwd = 1, range = 0, main = "", whisklty = 0, staplelty = 0)

                  }
                dev.off()

                png(file = ptpng, width = 800, height=800)
                nf <- layout(matrix(1:6,3,2,byrow = FALSE),
                             c(2.1,2), c(2,2,2), FALSE)
                par(mar=c(2,3.5,3.5,1))
                boxplot(rc[,1]~facx, col = colours[6], range =  0, main ="", whisklty = 0, staplelty = 0)
                mtext(sN[1], side = 3,padj=-1,adj=1.1,cex=1)
                mtext("Row", side = 3,padj=-0.5,cex=0.8)
                mtext("log(red intensity)", side = 2, adj = 0.5,padj=-4,cex=0.7)
                
                boxplot(gc[,1]~facx, col = colours[4], lwd = 1, range =  0, main = "", whisklty = 0, staplelty = 0)
                mtext("log(green intensity)", side = 2, adj = 0.5,padj=-4,cex=0.7)
                boxplot(dat[,1]~facx, col = colours[2], lwd = 1, range = 0, main = "", whisklty = 0, staplelty = 0)
                mtext("log(ratio)", side = 2, adj = 0.5,padj=-4,cex=0.7)
                
                boxplot(rc[,1]~facy, col = colours[6], range =  0, main ="", whisklty = 0, staplelty = 0)
                mtext("Column", side = 3,padj=-0.5,cex=0.8)
                
                boxplot(gc[,1]~facy, col = colours[4], lwd = 1, range =  0, main = "", whisklty = 0, staplelty = 0)
                boxplot(dat[,1]~facy, col = colours[2], lwd = 1, range = 0, main = "", whisklty = 0, staplelty = 0)                              
                dev.off()
                 
                pttext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S1.4\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", "Row-Column effect", basename(ptpdf), basename(ptpng), figure)

                legendpt = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the boxplots of the log<sub>2</sub> intensities grouped by row (left panel) and column (right panel) of the array. From the top to the bottom, red channel, green channel and log(ratio) are represented. If there is no spatial effect, the boxes should be homogeneous in size (IQR) and y position (median).</DIV>",  figure)
              }

############################################
###Section 2 : Homogeneity between arrays###
############################################

#############################
###Section 2.1 : Boxplots ###
#############################
            
            section = section + 1
            sec2text = sprintf("<hr><h2><a name = \"S2\">Section %s: Homogeneity between arrays</a></h2>", section)
            
            figure = figure + 1
            
            xlim = range(dat, na.rm=TRUE)
            xlimr = range(rc, na.rm=TRUE)
            xlimg = range(gc, na.rm=TRUE)
            ylimrg = c(min(c(xlimr[1],xlimg[1])), max(c(xlimr[2], xlimg[2])))

            colours = brewer.pal(12, "Paired")
            if(numArrays <= 50)
              {
                xname = sN
                mai = c(0,0.8,0.2,0.2)
                xaxt = "s"
              }
            if(numArrays > 50)
              {
                xname = FALSE
                mai = c(0,0.8,0.2,0.2)              
                xaxt = "n"
              }

            if(all(intgroup %in% names(phenoData(expressionset)@data)) && grouprep)
              {
                covar = pData(expressionset)[colnames(pData(expressionset))==intgroup[1]][,1]
                lev = levels(as.factor(covar))
                colourCovd = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[1])

                
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
                    boxplot(lredc, col = colourCovd[as.factor(covar)], las = 3, range = 0,
                            names = xname, ylim = ylimrg, ylab = "log(red intensity)", cex.lab = 1.2)
                    boxplot(lgreenc, col = colourCovd[as.factor(covar)], las = 3, range = 0,
                            names = xname, ylim = ylimrg, ylab = "log(green intensity)", cex.lab = 1.2)
                    boxplot(ldat, col = colourCovd[as.factor(covar)], las = 3, range = 0,
                            names = xname, ylim = xlim, ylab = "log(ratio)", cex.lab = 1.2)
                  }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</center></b></td></table>\n", title="Boxplots", fig = figure)
              } else {
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
                    boxplot(lredc, col = colours[6], las = 3, range = 0,
                            names = xname, ylim = ylimrg, ylab = "log(red intensity)", cex.lab = 1.2)
                    boxplot(lgreenc, col = colours[4], las = 3, range = 0,
                            names = xname, ylim = ylimrg, ylab = "log(green intensity)", cex.lab = 1.2)
                    boxplot(ldat, col = colours[2], las = 3, range = 0,
                            names = xname, ylim = xlim, ylab = "log(ratio)", cex.lab = 1.2)
                  }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</center></b></td></table>\n", title="Boxplots", fig = figure)
              }
            htmltext2 = mplot2[[2]]
            legendhom1 = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> presents boxplots of the log<sub>2</sub>(Intensities). Each box corresponds to one array. The left panel corresponds to the red channel. The middle panel shows the green channel. The right panel shows the boxplots of log<sub>2</sub>(ratio). Typically, one expects the boxes to have similar size (IQR) and y position (median).</DIV>", figure)          

      
############################
###Section 2.2 : Density ###
############################
            
            figure = figure + 1

            if(all(intgroup %in% names(phenoData(expressionset)@data)) && grouprep)
              {
                covar = pData(expressionset)[colnames(pData(expressionset))==intgroup[1]][,1]
                lev = levels(as.factor(covar))
                colourCovd = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[1])

                mplot3 = makePlot(con=con, name = "density",
                  w=10, h=10, fun = function() {
                    for(n in seq_len(max(group)))
                      {

                        nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),
                                     c(2.3,2,2), c(2,2), TRUE)
                        par(mar=c(0,5,1,1),xaxt = "n")
                        multi("density",lredc[group==n],xlimr,"Density","","", col=colourCovd[as.factor(covar[group==n])])
                        legend("topright",legend=lev,lty=1,col=colourCovd[as.factor(lev)], bty = "n", cex = 0.9)
                        par(mar=c(0,2,1,1),xaxt = "n")
                        multi("density",lgreenc[group==n],xlimg,"","","", col=colourCovd[as.factor(covar[group==n])])
                        legend("topright",legend=lev,lty=1,col=colourCovd[as.factor(lev)], bty = "n", cex = 0.9)
                        if(max(group) !=1)
                          mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -3)
                        par(mar=c(0,2,1,1),xaxt = "n")
                        multi("density",ldat[group==n],xlim,"","","", col=colourCovd[as.factor(covar[group==n])])
                        legend("topright",legend=lev,lty=1,col=colourCovd[as.factor(lev)], bty = "n", cex = 0.9)
                        par(mar=c(1,5,0,1), xaxt = "s")
                        multi("ecdf",lredc[group==n],xlimr,"ECDF","log(red intensity)","", col=colourCovd[as.factor(covar[group==n])])
                        par(mar=c(1,2,0,1), xaxt = "s")
                        multi("ecdf",lgreenc[group==n],xlimg,"","log(green intensity)","", col=colourCovd[as.factor(covar[group==n])])
                        par(mar=c(1,2,0,1), xaxt = "s")
                        multi("ecdf",ldat[group==n],xlim,"","log(ratio)","", col=colourCovd[as.factor(covar[group==n])])

                      }
                  }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", title="Density plots", fig = figure)
              } else {
                mplot3 = makePlot(con=con, name = "density",
                  w=10, h=10, fun = function() {
                    for(n in seq_len(max(group)))
                      {

                        nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),
                                     c(2.3,2,2), c(2,2), TRUE)
                        par(mar=c(0,5,1,1),xaxt = "n")
                        multi("density",lredc[group==n],xlimr,"Density","","")
                        legend("topright",legend= sN[group==n],lty=1,col=brewer.pal(9, "Set1"), bty = "n", cex = 1.1)
                        par(mar=c(0,2,1,1),xaxt = "n")
                        multi("density",lgreenc[group==n],xlimg,"","","")
                        legend("topright", legend=sN[group==n],lty=1,col=brewer.pal(9, "Set1"), bty = "n", cex = 1.1)
                        if(max(group) !=1)
                          mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -3)
                        par(mar=c(0,2,1,1),xaxt = "n")
                        multi("density",ldat[group==n],xlim,"","","")
                        legend("topright", legend=sN[group==n],lty=1,col=brewer.pal(9, "Set1"), bty = "n", cex = 1.1)
                        par(mar=c(1,5,0,1), xaxt = "s")
                        multi("ecdf",lredc[group==n],xlimr,"ECDF","log(red intensity)","")
                        par(mar=c(1,2,0,1), xaxt = "s")
                        multi("ecdf",lgreenc[group==n],xlimg,"","log(green intensity)","")
                        par(mar=c(1,2,0,1), xaxt = "s")
                        multi("ecdf",ldat[group==n],xlim,"","log(ratio)","")
                      }
                  }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", title="Density plots", fig = figure)
              }
            
            htmltext3 = mplot3[[2]]          
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
                ngc = c(2:9)
                colb = brewer.pal(9, "Blues")
                colg = brewer.pal(9, "Greens")
                colr = brewer.pal(9, "OrRd")         
                fac = round(as.numeric(as.matrix(featureData(expressionset)$GC)),-1)      
                gcpng = "GCcontent.png"
                gcpdf = "GCcontent.pdf"

                png(file = gcpng)
                nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9),3,3,byrow = TRUE),
                             widths =c(2.5,2,2), heights=c(1.9,1.9,2.9), TRUE)
                par(mar=c(0,5,1,1), xaxt = "n")
                multi("density",rc~fac,xlimr,"Density","","", col = colr[ngc])
                par(mar=c(0,2,1,1))
                multi("density",gc~fac,xlimg,"","","", col = colg[ngc])
                par(mar=c(0,2,1,1))
                multi("density",dat~fac,xlim,"","","", col = colb[ngc])

                par(mar=c(1,5,0,1), xaxt="s")
                multi("ecdf",rc~fac,xlimr,"ECDF","","", col = colr[ngc])
                par(mar=c(1,2,0,1))
                multi("ecdf",gc~fac,xlimg,"","","", col = colg[ngc])
                par(mar=c(1,2,0,1))
                multi("ecdf",dat~fac,xlim,"","","", col = colb[ngc])


                par(mar=c(5,5,2,1))
                boxplot(rc~fac, col = colr[ngc], range =  0, main ="", xlab="log(red intensity)")
                mtext("Boxplot",side = 2,adj = 0.5, padj = -4 ,cex = 0.8)
                par(mar=c(5,2,2,1))
                boxplot(gc~fac, col = colg[ngc], lwd = 1, range =  0, main = "", xlab="log(green intensity)")           
                par(mar=c(5,2,2,1))
                boxplot(dat~fac, col = colb[ngc], lwd = 1, range = 0, main = "", xlab="log(ratio)")
                dev.copy(pdf, file = gcpdf)
                dev.off()
                dev.off()

                gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b><td><center><a name = \"S3.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</b></center></tr></td></table>\n", "GC content effect  ", basename(gcpdf), basename(gcpng), figure)

                legendgc = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the distributions of the log<sub>2</sub> intensities grouped by the percentage of cytosines (C) and guanines (G) among the nucleotides in each probe. From the top to the bottom, kernel density estimates, empirical cumulative distribution functions (ECDF) and boxplots are represented. Box and line colors in the three panels correspond to the same groups, the darker is the colour, the higher is the GC content. Cytosine and guanine are able to form three hydrogen bonds, while adenine (A) and thymine (T) only form two, hence oligonucleotides with a higher proportion of C and G can form more stable hybridization bindings. This should result in higher intensities measured on the array, regardless of the abundance of target molecules.</DIV>",  figure)
              }

######################################
###Section 3.1 : Mapping of probes ###
######################################
    
            if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
              {
                pmap = probesmap(expressionset, numArrays, section, figure, dat, sN, xlim, type = 2, con = con)

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
            hmap = hmap(expressionset, sN, section, figure, outM, numArrays, intgroup)
            
            section = hmap$section
            figure = hmap$figure
            htmltext4 = hmap$htmltext4
            sec4text = hmap$sec4text
            legendheatmap = hmap$legendheatmap
            
##########################################
###Section 5 : Variance Mean Dependence###
##########################################
            
            ##meanSdplot
            msdp = msdp(expressionset, section, figure, con, dat)
            
            section = msdp$section
            figure = msdp$figure
            htmltext5 = msdp$htmltext5
            sec5text = msdp$sec5text
            legendsdmean = msdp$legendsdmean

##########################
### Scores computation ###
##########################

            sc = scores(expressionset=expressionset,numArrays=numArrays, M=M, ldat=ldat, outM=outM, dat=dat, rc=rc, gc=gc, maxc=maxc, maxr=maxr, nuse=NULL, rle=NULL)
            
            scores = matrix("",ncol=length(sc),nrow=numArrays)
            for( i in seq_len(length(sc)))
              scores[unlist(sc[[i]][seq_len(length(sc[[i]]))]),i]="*"
            
            
##########################
### Writing the report ###
##########################
            
            arg = as.list(match.call(expand.dots = TRUE))
            con = report(expressionset=expressionset, arg=arg, sNt=sNt, sN=sN, sec1text=sec1text, mapdf=mapdf, matext1=matext1, nfig=nfig, legendMA=legendMA, batext=batext, nfig2=nfig2, bapng=bapng, ftext=ftext, pttext=pttext, legendpt=legendpt, nfig3=nfig3, fpng=fpng, legendlocal=legendlocal, sec2text=sec2text, htmltext2=htmltext2, legendhom1=legendhom1, group=group, htmltext3=htmltext3, legendhom2=legendhom2, sec3text=sec3text, gctext=gctext, legendgc=legendgc, gotext=gotext, legendgo=legendgo, sec4text=sec4text, htmltext4=htmltext4, legendheatmap=legendheatmap, sec5text=sec5text, htmltext5=htmltext5, legendsdmean=legendsdmean, scores=scores)

            writeLines("</table>", con)
            z=sessionInfo("arrayQualityMetrics")
            version = z$otherPkgs[[1]]$Version
            rversion = sessionInfo()$R.version$version.string
            writeLines(sprintf("<hr><DIV style=\"font-size: 13; font-family: Lucida Grande\">This report has been created with arrayQualityMetrics %s under %s</DIV>",version, rversion), con)

            closeHtmlPage(con)
            
          })####end set method NChannelSet
          
         
#################################################################
#################################################################
#################### ExpressionSet Functions ####################
#################################################################
#################################################################

aqm.expressionset = function(expressionset, outdir = getwd(), force =
  FALSE, do.logtransform = FALSE, split.plots = FALSE,
    intgroup = "Covariate", arg, grouprep) 
  {
    olddir = getwd()
    on.exit(setwd(olddir))

    ##data preparation
    if(do.logtransform)
      {
        dat = logtransform(exprs(expressionset))
      } else dat = exprs(expressionset)
    
    dircreation(outdir, force)
   
    sN = sampleNames(expressionset)
    ##list to use multidensity, multiecdf and boxplot
    ldat = mat2list(dat)
    gN = featureNames(expressionset)
    numArrays = ncol(dat)
    
    ##second part of data preparation
    dprep = prepdata(sN, dat, numArrays, split.plots)
            
    long = dprep$long
    sN = dprep$sN
    sNt = dprep$sNt
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
    medArray = rowMedians(dat, na.rm=TRUE)
    M =  dat - medArray
    A =  (dat + medArray)/2

    MAplot = maplot(M, A, sN, numArrays)

    section = as.numeric(MAplot$section)
    figure = as.numeric(MAplot$figure)
    matext1 = as.character(MAplot$matext1)
    sec1text = as.character(MAplot$sec1text)
    nfig = as.numeric(MAplot$nfig)
    mapdf = as.character(MAplot$mapdf)   
    
    legendMA = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the MA plot for each array. M and A are defined as:<br>M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>)),<br>where I<sub>1</sub> is the intensity of the array studied and I<sub>2</sub> is the intensity of a \"pseudo\"-array, which contains for each probe the median of that probe's values over all arrays. Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A. Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).</DIV>", figure)     

#########################################
###Section 1.2 : Spatial Plots aBatch ###
#########################################
  
    if(is(expressionset, "AffyBatch"))
      {
        figure = figure + 1
        maxc = ncol(expressionset)
        maxr = nrow(expressionset)
        nfig3 = if(maxc*maxr < 1000000) ceiling(numArrays/6) else numArrays
        colourRamp <- rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))

        fpng = paste("foreground", seq_len(nfig3), ".png", sep="")
        fignu = 1
        aR = if(maxr > maxc) maxr/maxc else maxc/maxr

        for(a in seq_len(numArrays))
          {
            rfi = rank(dat[,a])
            mrfi = matrix(rfi,ncol=maxc,nrow=maxr,byrow=T)                                        
            if(maxr > maxc)
              mrfi = t(mrfi)
                
            if(maxc*maxr < 1000000){
              if(a %in% seq(1,numArrays,by=6))
                {
                  png(fpng[fignu], width = 500, height = 500*aR/2)
                  nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = TRUE),
                               widths = c(1.5,1.5,1.5),
                               heights = c(1.5*aR,1.5*aR),
                               respect = TRUE)
                }
            } else png(fpng[fignu], width = 500, height = 500*aR)
            
            par(xaxt = "n", yaxt = "n",mar=c(1,2,2,1))
            image(mrfi, col = colourRamp)
            mtext(sN[a],side = 3, adj = 0.5, padj = -1 ,cex = 0.7)
                
            if(maxc*maxr >= 1000000)
              {
                dev.off()
                fignu = fignu +1
              } else {            
                if((a%%6==0) || (a == numArrays))
                  {
                    dev.off()
                    fignu = fignu +1
                  }
              }
          }
                                   
        m = matrix(pretty(mrfi,9),nrow=1,ncol=length(pretty(mrfi,9)))
        llpng = "localisationlegend.png" 
        png(file= llpng, width = 3*72, height = 7*72)
        nf <- layout(1, widths = 0.9, heights = 3, respect = TRUE)
        image(m,xaxt="n",yaxt="n",ylab="Rank(Intensity)", col = colourRamp, cex.lab = 0.8)
        axis(2, label= as.list(pretty(mrfi,9)),at=seq(0,1,by=(1/(length(pretty(mrfi,9))-1))), cex.axis = 0.7, padj = 1)
        dev.off()
           
        ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name=\"S1.2\"><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></a></A><center><b>Figure %s</b></center></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td><td>\n", "Spatial distribution of feature intensities", basename(fpng[1]), basename(fpng[1]), figure, basename(llpng))

        legendlocal = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s:</b> False color representations of the arrays' spatial distributions of feature intensities. The color scale is shown in the panel on the right, and it is proportional to the ranks of the probe intensities. This has indeed the potential to detect patterns that are small in amplitude, but systematic within array. These may not be consequential for the downstream data analysis, but if properly interpreted, could e.g. still be useful for technology and experimental protocol optimisation as it helps in identifying patterns that may be caused by, for example, spatial gradients in the hybridization chamber, air bubbles, spotting or plating problems.</DIV>", figure)

######################################
###Section 1.4 : R/C effect aBatch ###
######################################
        
        if(maxc*maxr < 1000000)
          {
            figure = figure + 1
            ptpdf = "row_column.pdf"
            ptpng = "row_column.png"
            colours = brewer.pal(12, "Paired")       
       
            pdf(file = ptpdf, width=10, height=4)
            nf <- layout(matrix(1:2,1,2,byrow = FALSE),
                         c(2,2), 2, FALSE)
        
            for ( a in seq_len(numArrays))
              {
                matdat = matrix(dat[,a],ncol=maxc,nrow=maxr,byrow=T)
                par(mar=c(2,3.5,3.5,1), cex.axis = 0.7)
                boxplot(matdat~row(matdat), col = colours[2], range =  0, main ="", whisklty = 0, staplelty = 0)
                mtext(sN[a], side = 3,padj=-1,adj=1.1,cex=1)
                mtext("Row", side = 3,padj=-0.5,cex=0.8)
                boxplot(matdat~col(matdat), col = colours[2], range =  0, main ="", whisklty = 0, staplelty = 0)
                mtext("Column", side = 3,padj=-0.5,cex=0.8)
              }
            dev.off()
        
            png(file = ptpng, width = 1000, height=400)
            nf <- layout(matrix(1:2,1,2,byrow = FALSE),
                         c(2,2), 2, FALSE)
            matdat = matrix(dat[,1],ncol=maxc,nrow=maxr,byrow=T)
            par(mar=c(2,3.5,3.5,1), cex.axis = 0.7)
            boxplot(matdat~row(matdat), col = colours[2], range =  0, main ="", whisklty = 0, staplelty = 0)
            mtext(sN[a], side = 3,padj=-1,adj=1.1,cex=1)
            mtext("Row", side = 3,padj=-0.5,cex=0.8)
            boxplot(matdat~col(matdat), col = colours[2], range =  0, main ="", whisklty = 0, staplelty = 0)
            mtext("Column", side = 3,padj=-0.5,cex=0.8)
            dev.off()
        
            pttext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S1.4\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", "Row-Column effect", basename(ptpdf), basename(ptpng), figure)
        
            legendpt = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the boxplots of the log<sub>2</sub> intensities grouped by row (left panel) and column (right panel) of the array. Typically, one expects the boxes to have similar size (IQR) and y position (median).</DIV>",  figure)
          }      
      }

#######################################
###Section 1.2 : Spatial Plots Eset ###
#######################################
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
            fpng = paste("foreground", seq_len(nfig3), ".png", sep="")
            fignu = 1
            b = 1
    
            aR = if(maxr > maxc) maxr/maxc else maxc/maxr

            intf = cbind(r,c,dat)
            for(a in seq_len(numArrays))
              {
                fg = matrix(NA,ncol=maxc,nrow=maxr)
                   
                for(i in seq_len(nrow(intf)))
                  fg[intf[i,1],intf[i,2]] = intf[i,(2+a)]
                    
                mfg = matrix(rank(fg),ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
                    
                if(maxr > maxc)
                  mfg = t(mfg)
                

                if(a %in% seq(1,numArrays,by=6))
                  {
                    png(fpng[fignu], width = 350, height = 350*aR/2)
                    nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = FALSE),
                                 respect = FALSE)
                  }

                par(xaxt = "n", yaxt = "n",mar=c(1,2,2,1))
                image(mfg, col = colourRamp)
                mtext(sN[a],side = 3, line = 0.5 ,cex = 0.7)
                if((a%%6 == 0) || (a == numArrays))                   
                  {
                    dev.off()
                    fignu = fignu +1
                  }
              }
                              
            m = matrix(pretty(mfg,9),nrow=1,ncol=length(pretty(mfg,9)))
            llfpng = "localisationlegendforeground.png"
            png(file= llfpng, width = 3*72, height = 7*72)
            nf <- layout(1, widths = 0.9, heights = 3, respect = TRUE)

            image(m,xaxt="n",yaxt="n",ylab="Rank(Intensity)", col = colourRamp, cex.lab = 0.8)
            axis(2, label= as.list(pretty(mfg,9)),at=seq(0,1,by=(1/(length(pretty(mfg,9))-1))), cex.axis = 0.7, padj = 1)
            dev.off()
                
            ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name = \"S1.3\"><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></a></A><center><b>Figure %s</b></center></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td><td>\n", "Spatial distribution of feature intensites", basename(fpng[1]), basename(fpng[1]), figure, basename(llfpng))

            legendlocal = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s:</b> False color representations of the arrays' spatial distributions of feature intensities. The color scale is shown in the panel on the right, and it is proportional to the ranks of the probe intensities. This has indeed the potential to detect patterns that are small in amplitude, but systematic within array. These may not be consequential for the downstream data analysis, but if properly interpreted, could e.g. still be useful for technology and experimental protocol optimisation as it helps in identifying patterns that may be caused by, for example, spatial gradients in the hybridization chamber, air bubbles, spotting or plating problems.</DIV>", figure)

#####################################
###Section 1.4 : R/C effect  Eset ###
#####################################
            
            figure = figure + 1
            ptpdf = "rowcolumn.pdf"
            ptpng = "rowcolumn.png"
            colours = brewer.pal(12, "Paired")
            
            facx = round(as.numeric(as.matrix(featureData(expressionset)$X)),-1)
            facy = round(as.numeric(as.matrix(featureData(expressionset)$Y)),-1)         
            pdf(file = ptpdf, width=10, height=4)
            nf <- layout(matrix(1:2,1,2,byrow = FALSE),
                         c(2,2), 2, FALSE)
            
            for ( a in seq_len(numArrays))
              {
                par(mar=c(2,3.5,3.5,1), cex.axis = 0.7)
                boxplot(dat[,a]~facx, col = colours[2], range =  0, main ="", whisklty = 0, staplelty = 0)
                mtext(sN[a], side = 3,padj=-1,adj=1.1,cex=1)
                mtext("Row", side = 3,padj=-0.5,cex=0.8)
                boxplot(dat[,a]~facy, col = colours[2], range =  0, main ="", whisklty = 0, staplelty = 0)
                mtext("Column", side = 3,padj=-0.5,cex=0.8)
              }
            dev.off()
            
            png(file = ptpng, width = 1000, height=400)
            nf <- layout(matrix(1:2,1,2,byrow = FALSE),
                         c(2,2), 2, FALSE)
            par(mar=c(2,3.5,3.5,1), cex.axis=0.7)
            boxplot(dat[,1]~facx, col = colours[2], range =  0, main ="", whisklty = 0, staplelty = 0)
            mtext(sN[1], side = 3,padj=-1,adj=1.1,cex=1)
            mtext("Row", side = 3,padj=-0.5,cex=0.8)
            boxplot(dat[,1]~facy, col = colours[2], range =  0, main ="", whisklty = 0, staplelty = 0)
            mtext("Column", side = 3,padj=-0.5,cex=0.8)
            dev.off()
            
            pttext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S1.4\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", "Row-Column effect", basename(ptpdf), basename(ptpng), figure)

            legendpt = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the boxplots of the log<sub>2</sub> intensities grouped by row (left panel) and column (right panel) of the array.  Typically, one expects the boxes to have similar size (IQR) and y position (median).</DIV>",  figure)
          }
      }

############################################
###Section 2 : Homogeneity between arrays###
############################################

#############################
###Section 2.1 : Boxplots ###
#############################

    section = section + 1
    sec2text = sprintf("<hr><h2><a name = \"S2\">Section %s: Homogeneity between arrays</h2></a>", section)
            
    figure = figure + 1
    if(numArrays <= 50)
      {
        xname = sN
        xaxt = "s"
      }
    if(numArrays > 50)
      {
        xname = FALSE
        xaxt = "n"
      }
    if(all(intgroup %in% names(phenoData(expressionset)@data)) && grouprep)
      {
        covar = pData(expressionset)[colnames(pData(expressionset))==intgroup[1]][,1]
        lev = levels(as.factor(covar))
        colourCovd = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[1])

        mplot2 = makePlot(con=con, name = "boxplot",
          w=8, h=8, fun = function() {
            par(cex.axis = 1, pty = "s", lheight =((1/log10(numArrays))*long), mai = c(((long/12)+0.2),0.4,0.2,0.2) , omi = c(0,0,0,0))
            boxplot(ldat, col = colourCovd[as.factor(covar)], las = 3, range = 0, names = xname)
          }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</center></b></td></table>\n", title = "Boxplots", fig = figure)
      } else {
        mplot2 = makePlot(con=con, name = "boxplot",
          w=8, h=8, fun = function() {
            colours = brewer.pal(12, "Paired")
            par(cex.axis = 1, pty = "s", lheight =((1/log10(numArrays))*long), mai = c(((long/12)+0.2),0.4,0.2,0.2) , omi = c(0,0,0,0))
            boxplot(ldat, col = colours[2], las = 3, range = 0, names = xname)
          }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</center></b></td></table>\n", title = "Boxplots", fig = figure)
      }
    
    htmltext2 = mplot2[[2]]
    legendhom1 = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> presents boxplots of the log<sub>2</sub>(Intensities). Each box corresponds to one array. Typically, one expects the boxes to have similar size (IQR) and y position (median).</DIV>", figure)          
    
############################
###Section 2.2 : Density ###
############################
    
    figure = figure + 1
    xlim = range(dat, na.rm=TRUE)

    if(all(intgroup %in% names(phenoData(expressionset)@data)) && grouprep)
      {
        covar = pData(expressionset)[colnames(pData(expressionset))==intgroup[1]][,1]
        lev = levels(as.factor(covar))
        colourCovd = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[1])

        mplot3 = makePlot(con=con, name = "density",
          w=10, h=10, fun = function() {
            for(n in seq_len(max(group)))
              {            
                nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(2.8,2.8),c(1.8,2), TRUE)
                par(xaxt = "n", mar = c(0,5,2,5))
                multi("density",ldat[group==n],xlim,"Density","","", col=colourCovd[as.factor(covar[group==n])])
                legend("topright",legend=lev,lty=1,col=colourCovd[as.factor(lev)], bty = "n", cex = 0.9)
                if(max(group) != 1)
                  mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -1)
                par(xaxt = "s", mar = c(4,5,0,5))
                multi("ecdf",ldat[group==n],xlim,"ECDF","log(intensity)","", col=colourCovd[as.factor(covar[group==n])])
              }} , text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", title="Density plots", fig = figure)
      } else{    
        mplot3 = makePlot(con=con, name = "density",
          w=10, h=10, fun = function() {
            for(n in seq_len(max(group)))
              {            
                nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(2.8,2.8),c(1.8,2), TRUE)
                par(xaxt = "n", mar = c(0,5,2,5))
                multi("density",ldat[group==n],xlim,"Density","","")
                legend("topright", legend=sN[group==n],lty=1,col=brewer.pal(9, "Set1"), bty = "n", cex = 0.9)
                if(max(group) != 1)
                  mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -1)
                par(xaxt = "s", mar = c(4,5,0,5))
                multi("ecdf",ldat[group==n],xlim,"ECDF","log(intensity)","")
              }} , text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", title="Density plots", fig = figure)
      }
    htmltext3 = mplot3[[2]]
    
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
        ngc = c(2:9)
        colb = brewer.pal(9, "Blues")
        colg = brewer.pal(9, "Greens")
        colr = brewer.pal(9, "OrRd")         
        fac = round(as.numeric(as.matrix(featureData(expressionset)$GC)),-1)
               
        gcpng = "GCcontent.png"
        gcpdf = "GCcontent.pdf"
        png(file = gcpng)
        nf <- layout(matrix(c(1,2,3),3,1,byrow = FALSE),
                     2.2, c(2,1.8,2), TRUE)
        par(mar=c(0,3.5,3.5,1))
        multi("density",dat~fac,xlim,"Density","","Log(ratio)", col = colb[ngc],  xaxt = "n")
        par(mar=c(2,3.5,0,1))
        multi("ecdf",dat~fac,xlim,"ECDF","","", col = colb[ngc])
        par(mar=c(2,3.5,2,1))
        boxplot(dat~fac, col = colb[ngc], lwd = 1, range = 0)
        mtext("Boxplot",side = 2,adj = 0.5, padj = -3.3 ,cex = 0.7)
        dev.copy(pdf, file = gcpdf)
        dev.off()
        dev.off()
                
        gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b><td><center><a name = \"S3.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</b></center></tr></td></table>\n", "GC content effect  ", basename(gcpdf), basename(gcpng), figure)

        legendgc = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the distributions of the log<sub>2</sub> intensities grouped by the percentage of cytosines (C) and guanines (G) among the nucleotides in each probe. From the top to the bottom, kernel density estimates, empirical cumulative distribution functions (ECDF) and boxplots are represented. Box and line colors in the three panels correspond to the same groups. Cytosine and guanine are able to form three hydrogen bonds, while adenine (A) and thymine (T) only form two, hence oligonucleotides with a higher proportion of C and G can form more stable hybridization bindings. This should result in higher intensities measured on the array, regardless of the abundance of target molecules.</DIV>",  figure)
      }

######################################
###Section 3.2 : Mapping of probes ###
######################################
    
    if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
      {
        pmap = probesmap(expressionset, numArrays, section, figure, dat, sN, xlim, type = 1, con = con)

        section = pmap$section
        figure = pmap$figure
        gotext = pmap$gotext
        sec3text = as.character(pmap$sec3text)
                
        legendgo = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the density distributions of the log<sub>2</sub> intensities grouped by the mapping of the probes. Blue, density estimate of intensities of probes annotated \"TRUE\" in the <b>\"hasTarget\"</b> slot. Gray, probes annotated \"FALSE\" in the <b>\"hasTarget\"</b> slot. We expect that, probes mapping with a coding mRNA will have higher expression level than unannotated probes, the density curve of mapping probes should be shifted to the right of the curve of unannotated probes.</DIV>", figure)          
      }
    
##########################################
###Section 4 : Between array comparison###
##########################################
   
    ##Heatmap
    hmap = hmap(expressionset, sN, section, figure, outM, numArrays, intgroup)            
    section = hmap$section
    figure = hmap$figure
    htmltext4 = hmap$htmltext4
    sec4text = hmap$sec4text
    legendheatmap = hmap$legendheatmap

##########################################
###Section 5 : Variance Mean Dependence###
##########################################
    
    ##meanSdplot
    msdp = msdp(expressionset, section, figure, con, dat)            
    section = msdp$section
    figure = msdp$figure
    htmltext5 = msdp$htmltext5
    sec5text = msdp$sec5text
    legendsdmean = msdp$legendsdmean

    
##########################
### Scores computation ###
##########################
    if(is(expressionset, "ExpressionSet"))
      {
        sc = scores(expressionset=expressionset,numArrays=numArrays,
          M=M, ldat=ldat, outM=outM, dat=dat, rc=NULL, gc=NULL,
          maxc=maxc, maxr=maxr, nuse=NULL, rle=NULL)
        
        scores = matrix("",ncol=length(sc),nrow=numArrays)
        for( i in seq_len(length(sc)))
          scores[unlist(sc[[i]][seq_len(length(sc[[i]]))]),i]="*"
      }


#########################
### Writing the report###
#########################
    
    if(is(expressionset, "ExpressionSet"))
      {
        con = report(expressionset=expressionset, arg=arg, sNt=sNt,
          sN=sN, sec1text=sec1text, mapdf=mapdf, matext1=matext1,
          nfig=nfig, legendMA=legendMA,ftext=ftext, pttext=pttext,
          legendpt=legendpt, nfig3=nfig3, fpng=fpng,
          legendlocal=legendlocal, sec2text=sec2text, htmltext2=htmltext2,
          legendhom1=legendhom1, group=group, htmltext3=htmltext3,
          legendhom2=legendhom2, sec3text=sec3text, gctext=gctext,
          legendgc=legendgc, gotext=gotext, legendgo=legendgo,
          sec4text=sec4text, htmltext4=htmltext4,
          legendheatmap=legendheatmap, sec5text=sec5text,
          htmltext5=htmltext5, legendsdmean=legendsdmean, scores=scores)
        
        l = list(numArrays=numArrays,sN=sN, section=section,
          figure=figure, con=con, dat=dat, olddir=olddir) 
      }
    
    if(is(expressionset, "AffyBatch"))
      {
        l = list(numArrays=numArrays, section=section, figure=figure,
          dat=dat, olddir=olddir, expressionset=expressionset, arg=arg,
          sN=sN, sec1text=sec1text, mapdf=mapdf, matext1=matext1,
          nfig=nfig, legendMA=legendMA, ftext=ftext, nfig3=nfig3,
          fpng=fpng, legendlocal=legendlocal, sec2text=sec2text,
          htmltext2=htmltext2, legendhom1=legendhom1, group=group,
          htmltext3=htmltext3, legendhom2=legendhom2,  sec4text=sec4text,
          htmltext4=htmltext4, legendheatmap=legendheatmap,
          sec5text=sec5text, htmltext5=htmltext5,
          legendsdmean=legendsdmean, M=M, outM=outM, maxc=maxc, maxr=maxr,
          ldat=ldat) 
        l$sNt = sNt
        if("GC" %in% rownames(featureData(expressionset)@varMetadata))
          {
            l$sec3text=sec3text
            l$gctext=gctext
            l$legendgc=legendgc
          }
        if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
          {
            l$sec3text=sec3text
            l$gotext=gotext
            l$legendgo=legendgo
          }       
        if(maxc*maxr < 1000000)
          {
            l$pttext=pttext
            l$legendpt=legendpt
          }       
      }
    
    return(l) 
}

#################################################################
#################################################################
###################### ExpressionSet Method #####################
#################################################################
#################################################################

setMethod("arrayQualityMetrics",signature(expressionset="ExpressionSet"),
          function(expressionset, outdir, force, do.logtransform,
          split.plots, intgroup, grouprep) 
          {
            arg = as.list(match.call(expand.dots = TRUE))

            l = aqm.expressionset(expressionset, outdir, force,
                do.logtransform, split.plots, intgroup, arg, grouprep) 
            con = l$con
            olddir = l$olddir
            on.exit(setwd(olddir))
            setwd(outdir)

            
            writeLines("</table>", con)
            z=sessionInfo("arrayQualityMetrics")
            version = z$otherPkgs[[1]]$Version
            rversion = sessionInfo()$R.version$version.string
            writeLines(sprintf("<hr><DIV style=\"font-size: 13; font-family: Lucida Grande\">This report has been created with arrayQualityMetrics %s under %s</DIV>",version, rversion), con)
            
            closeHtmlPage(con)
          }) ##end set method ExpressionSet

#################################################################
#################################################################
######################## AffyBatch Method #######################
#################################################################
#################################################################

setMethod("arrayQualityMetrics",signature(expressionset="AffyBatch"),
          function(expressionset, outdir, force, do.logtransform,
          split.plots, intgroup, grouprep)
          {
            
            arg = as.list(match.call(expand.dots = TRUE))

            l = aqm.expressionset(expressionset, outdir, force, do.logtransform, split.plots, intgroup, arg, grouprep)
            numArrays = as.numeric(l$numArrays)
            expressionset = l$expressionset
            sN = l$sN
            sNt = l$sNt
            section = as.numeric(l$section)
            figure = as.numeric(l$figure)
            dat = l$dat
            olddir = l$olddir
            setwd(outdir)

            on.exit(setwd(olddir))


############################
###Section 8 : Affy plots###
############################

            cols = brewer.pal(9, "Set1")            
            section = section + 1
            sec6text = sprintf("<hr><h2><a name = \"S6\">Section %s: Affymetrix specific plots</h2></a>", section)
            
            figure1 = figure + 1
            acol = sample(brewer.pal(8, "Dark2"), numArrays, replace = (8<numArrays))
            rnaDeg = try(AffyRNAdeg(expressionset, log.it = do.logtransform))
            if(class(rnaDeg)=='try-error')
              warning("RNA degradation plot from the package 'affy' cannot be produced for this data set.")
            affypng1 = "RNAdeg.png"
            affypdf1 = "RNAdeg.pdf"
            png(file = affypng1)
            try(plotAffyRNAdeg(rnaDeg, cols = acol, lwd = 2))
            dev.copy(pdf, file = affypdf1)
            dev.off()
            dev.off()
                
            figure2 = figure1 + 1
            pp1 = try(preprocess(expressionset))
            dataPLM = try(fitPLM(pp1, background = FALSE, normalize = FALSE))
            if(class(pp1)=='try-error' || class(dataPLM)=='try-error')
              warning("RLE and NUSE plots from the package 'affyPLM' cannot be produced for this data set.")
            affypng2 = "RLE.png"
            affypdf2 = "RLE.pdf"
            png(file = affypng2)
            rle = try(Mbox(dataPLM, ylim = c(-1, 1),
              names = sN, col = cols[2],
              whisklty = 0, staplelty = 0,
              main = "RLE", las = 3, cex.axis = 0.8))
            dev.copy(pdf, file = affypdf2)
            dev.off()
            dev.off()
            
            figure3 = figure2 + 1
            affypng3 = "NUSE.png"
            affypdf3 = "NUSE.pdf"
            png(file = affypng3)
            nuse = try(boxplot(dataPLM, ylim = c(0.95, 1.5), names = sN,
              outline = FALSE, col = cols[2], main = "NUSE",
              las = 2, cex.axis = 0.8))
            dev.copy(pdf, file = affypdf3)
            dev.off()
            dev.off()
            
            qcStats = try(qc(expressionset))
            figure4 = figure3 + 1
            affypng4 = "qc.png"
            affypdf4 = "qc.pdf"
            png(file = affypng4)
            p = try(plot(qcStats, cex.axis = 0.8))
            suppressWarnings(dev.copy(pdf, file = affypdf4))
            dev.off()
            dev.off()
            if(class(p)=='try-error')
              warning("QCstat plot from the package 'simpleaffy' cannot be produced for this data set.")

            
            affytext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name = \"S6.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><center><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><a name = \"S6.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><center><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><a name = \"S6.3\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><center><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><a name = \"S6.4\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><center><BR><b>Figure %s</b></CENTER></tr></td></table>\n", "RNA degradation plot", basename(affypdf1), basename(affypng1), figure1, "RLE plot", basename(affypdf2), basename(affypng2), figure2, "NUSE plot", basename(affypdf3), basename(affypng3), figure3, "Diagnostic plot recommended by Affymetrix", basename(affypdf4), basename(affypng4), figure4)
            legendaffy = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\">In this section we present diagnostic plots based on tools provided in the affyPLM package. In <b>Figure %s</b> a RNA digestion plot is computed on normalized data (so that standard deviation is equal to 1). In this plot each array is represented by a single line. It is important to identify any array(s) that has a slope which is very different from the others. The indication is that the RNA used for that array has potentially been handled quite differently from the other arrays.  <b>Figure %s</b> is a Relative Log Expression (RLE) plot and an array that has problems will either have larger spread, or will not be centered at M = 0, or both. <b>Figure %s</b> is a Normalized Unscaled Standard Error (NUSE) plot. Low quality arrays are those that are substantially elevated or more spread out, relative to the other arrays. NUSE values are not comparable across data sets. Both RLE and NUSE are performed on preprocessed data (background correction and quantile normalization). <b>Figure %s</b> represents the diagnostic plot recommended by Affymetrix. It is fully describe in the simpleaffy.pdf vignette of the package simpleaffy. Any metrics that is shown in red is out of the manufacturer's specific boundaries and suggests a potential problem, any metrics shown in blue is fine.</DIV>", figure1, figure2, figure3, figure4)            
            
            ##PM.MM
            figure5 = figure4 + 1

            cols = brewer.pal(9, "Set1")
            xlim = range(dat, na.rm=TRUE)
            pmopng = "PM.MM.png"
            pmopdf = "PM.MM.pdf"       
        
            png(file = pmopng)
            pm = try(pmmm(expressionset,xlim,"","","", cex.axis = 0.9))
            if(class(pm)=='try-error')
              {
                warning("PM MM plot cannot be produced for this data set.")
                dev.off()
                pmotext = NULL
                legendpmo = NULL
              }

            if(class(pm)!='try-error')
              {
                legend("topright", c("PM","MM"),lty=1,lwd=2,col= c(cols[c(2,9)]), bty = "n")
                dev.copy(pdf, file = pmopdf)
                dev.off()
                dev.off()
                         
                pmotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b><td><center><a name = \"S6.5\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></A><br><b>Figure %s</b></center></tr></td></table>\n", "Perfect matchs and mismatchs ", basename(pmopdf), basename(pmopng), figure5)            
           
                legendpmo = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the density distributions of the log<sub>2</sub> intensities grouped by the matching of the probes. Blue, density estimate of intensities of perfect match probes (PM) and gray the mismatch probes (MM). We expect that, MM probes having poorer hybridization than PM probes, the PM curve should be shifted to the right of the MM curve.</DIV>",  figure5)
              }

##########################
### Scores computation ###
##########################
            sc = scores(expressionset=expressionset,numArrays=numArrays, M=l$M, ldat=l$ldat, outM=l$outM, dat=dat, rc=NULL, gc=NULL, maxc=l$maxc, maxr=l$maxr, nuse=nuse, rle=rle)               
            scores = matrix("",ncol=length(sc),nrow=numArrays)
            for( i in seq_len(length(sc)))
              scores[unlist(sc[[i]][seq_len(length(sc[[i]]))]),i]="*"
            
            
#########################
### Writing the report###
#########################
            con = report(expressionset=expressionset, arg=arg,
              sNt=sNt, sN=sN, sec1text=l$sec1text, mapdf=l$mapdf,
              matext1=l$matext1, nfig=l$nfig, legendMA=l$legendMA,
              batext=l$batext, nfig2=l$nfig2, bapng=l$bapng,
              ftext=l$ftext, pttext=l$pttext, legendpt=l$legendpt,
              nfig3=l$nfig3, fpng=l$fpng, legendlocal=l$legendlocal,
              sec2text=l$sec2text, htmltext2=l$htmltext2,
              legendhom1=l$legendhom1, group=l$group,
              htmltext3=l$htmltext3, legendhom2=l$legendhom2,
              sec3text=l$sec3text, gctext=l$gctext,
              legendgc=l$legendgc, gotext=l$gotext,
              legendgo=l$legendgo, sec4text=l$sec4text,
              htmltext4=l$htmltext4, legendheatmap=l$legendheatmap,
              sec5text=l$sec5text, htmltext5=l$htmltext5,
              legendsdmean=l$legendsdmean, scores=scores)
            
            writeLines(sec6text, con)
            writeLines(affytext, con)
            writeLines(legendaffy, con)
            if(!is.null(pmotext))
              {
                writeLines(pmotext, con)
                writeLines(legendpmo, con)
              }


            writeLines("</table>", con)
            z=sessionInfo("arrayQualityMetrics")
            version = z$otherPkgs[[1]]$Version
            rversion = sessionInfo()$R.version$version.string
            writeLines(sprintf("<hr><DIV style=\"font-size: 13; font-family: Lucida Grande\">This report has been created with arrayQualityMetrics %s under %s</DIV>",version, rversion), con)

            closeHtmlPage(con)

          }) ##end functions for AffyBatch and ExpressionSet

###################################################################
###################################################################
#######################  BeadLevelList Method #####################
###################################################################
###################################################################

setMethod("arrayQualityMetrics",signature(expressionset = "BeadLevelList"),
          function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            olddir = getwd()
            on.exit(setwd(olddir))

            dircreation(outdir, force)
            
            sN = arrayNames(expressionset)
            numArrays = as.numeric(dim(expressionset)[1])
            
            if(expressionset@arrayInfo$channels == "single")
              summaryES = createBeadSummaryData(expressionset, imagesPerArray = 1, log = do.logtransform)
            if(expressionset@arrayInfo$channels == "two")
              {
                summaryESRG = createBeadSummaryData(expressionset, what = "RG", imagesPerArray = 1, log = do.logtransform)
                rc = assayData(summaryESRG)$R                         
                gc = assayData(summaryESRG)$G                
                summaryES = createBeadSummaryData(expressionset ,what = "M", imagesPerArray = 1, log = do.logtransform)
              }            

            dat = exprs(summaryES)
            
            ##second part of data preparation
            dprep = prepdata(sN, dat, numArrays, split.plots)
            long = dprep$long
            sN = dprep$sN
            sNt = dprep$sNt
            outM = dprep$outM
            group = dprep$group

            ##one or two colour settings
            switch(expressionset@arrayInfo$channels,
                   "single" = {
                     ndiv = 6
                     dispo = TRUE
                     widthbox = 8
                     d = lapply(seq_len(numArrays), function(i) getArrayData(expressionset, array = i, what = "G",log = do.logtransform))
                     xlim = range(unlist(d), na.rm=TRUE)
                   },
                   "two" = {
                     ndiv = 3
                     dispo = FALSE
                     widthbox = 15
                     dr  = lapply(seq_len(numArrays), function(i) getArrayData(expressionset, array = i, what = "R",log = do.logtransform))
                     dg  = lapply(seq_len(numArrays), function(i) getArrayData(expressionset, array = i, what = "G",log = do.logtransform))
                     dlr = lapply(seq_len(numArrays), function(i) getArrayData(expressionset, array = i, what = "M",log = do.logtransform))
                     xlim  = range(unlist(dlr), na.rm=TRUE)
                     xlimr = range(unlist(dr), na.rm=TRUE)
                     xlimg = range(unlist(dg), na.rm=TRUE)
                   },
                   stop(sprintf("Invalid expressionset@arrayInfo$channels '%s'.", expressionset@arrayInfo$channels)))


####################################
###Section 1 : Per array quality ###
####################################       
           
#############################
###Section 1.1 : MA-plots ###
#############################            
            ##MA-plots
            if(expressionset@arrayInfo$channels == "single")
              {
                medArray = rowMedians(dat, na.rm=TRUE)
                M =  dat - medArray
                A =  (dat + medArray)/2
              }
            if(expressionset@arrayInfo$channels == "two")
              {
                M = rc - gc
                A = 0.5*(rc +gc)
              }
            
            MAplot = maplot(M, A, sN, numArrays)
            
            section = as.numeric(MAplot$section)
            figure = as.numeric(MAplot$figure)
            matext1 = as.character(MAplot$matext1)
            sec1text = as.character(MAplot$sec1text)
            nfig = as.numeric(MAplot$nfig)
            mapdf = as.character(MAplot$mapdf)

            if(expressionset@arrayInfo$channels == "single")
              legspe = "where I<sub>1</sub> is the intensity of the array studied and I<sub>2</sub> is the intensity of a \"pseudo\"-array, which contains for each probe the median of that probe's values over all arrays."
            if(expressionset@arrayInfo$channels == "two")
              legspe = "where I<sub>1</sub> and I<sub>2</sub> are the vectors of intensities of the two channels."
            ## fixme: consider using 'switch' in such situations
            
            legendMA = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> represents MA plot for each array. M and A are defined as :<br>
M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>
A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>))<br>
%s The calculations are done on the summarized data obtained by using the function createBeadSummaryData from the package beadarray. Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A. Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).</DIV>", figure, legspe)                   
            
#################################
###Section 1.2 : Spatial plots###
#################################

            figure = figure + 1            
            nfig2 = ceiling(numArrays/ndiv)      
            
            bapng = paste("background", seq_len(nfig2), ".png", sep="")
            fignu = 1
            
            for(a in seq_len(numArrays))
              {               
                if(a %in% seq(1,numArrays,by=ndiv))
                  {
                    png(bapng[fignu], width = 500, height = 500)
                    nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = dispo),
                                 respect = TRUE)
                  }
                
                if(expressionset@arrayInfo$channels == "single")
                  {
                    imageplot(expressionset, array = a, high = "yellow", low = "blue", main="", whatToPlot="Gb", log = do.logtransform)
                    mtext(sN[a],side = 3, adj = 0.5, padj = -1 ,cex = 0.7)
                  }
                
                if(expressionset@arrayInfo$channels == "two")
                  {
                    imageplot(expressionset, array = a, high = "yellow", low = "blue", main="", whatToPlot="Rb", log = do.logtransform)
                    mtext(sN[a],side = 3,adj = 0.5, padj = -1 ,cex = 0.7)
                    mtext("Red Intensity",side = 2, line = 0.1 ,cex = 0.7)
                    imageplot(expressionset, array = a, high = "yellow", low = "blue", main="", whatToPlot="Gb", log = do.logtransform)
                    mtext("Green Intensity",side = 2, line = 0.1 ,cex = 0.7)
                  }
                
                if((a%%ndiv == 0) || (a == numArrays))                   
                  {
                    dev.off()
                    fignu = fignu +1
                  }
              }
            
            batext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S1.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A></CENTER></td><td>\n", "Spatial distribution of local background intensites", basename(bapng[1]), basename(bapng[1]))


###################################
###Section 1.3 : Spatial plots 2###
###################################

            figure = figure + 1
            nfig3 = ceiling(numArrays/ndiv)      
            
            fpng = paste("foreground", seq_len(nfig3), ".png", sep="")
            fignu = 1
            
            for(a in seq_len(numArrays))
              {               
                if(a %in% seq(1,numArrays,by=ndiv))
                  {
                    png(fpng[fignu], width = 500, height = 500)
                    nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = dispo),
                                 respect = TRUE)
                  }
                
                if(expressionset@arrayInfo$channels == "single")
                  {
                    imageplot(expressionset, array = a, high = "yellow", low = "blue", main="", whatToPlot="G", log = do.logtransform)
                    mtext(sN[a],side = 3, adj = 0.5, padj = -1 ,cex = 0.7)
                  }

                if(expressionset@arrayInfo$channels == "two")
                  {
                    imageplot(expressionset, array = a, high = "yellow", low = "blue", main="", whatToPlot="R", log = do.logtransform)
                    mtext(sN[a],side = 3,adj = 0.5, padj = -1 ,cex = 0.7)
                    mtext("Red Intensity",side = 2, line = 0.1 ,cex = 0.7)
                    imageplot(expressionset, array = a, high = "yellow", low = "blue", main="", whatToPlot="G", log = do.logtransform)
                    mtext("Green Intensity",side = 2, line = 0.1 ,cex = 0.7)
                  }
               
                if((a%%ndiv == 0) || (a == numArrays))                   
                  {
                    dev.off()
                    fignu = fignu +1
                  }
              }
                        
            ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><a name=\"S1.3\"><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></a></A><center><b>Figure %s</b></center></td><td>\n", "Spatial distribution of feature intensities", basename(fpng[1]), basename(fpng[1]), figure)
            
            legendlocal = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s:</b> False color representations of the arrays' spatial distributions of feature intensities and, if available, local background estimates. These plots may help in identifying patterns that may be caused by, for example, spatial gradients in the hybridization chamber, air bubbles, spotting or plating problems.</DIV>", figure)          


############################################
###Section 2 : Homogeneity between arrays###
############################################

#############################
###Section 2.1 : Boxplots ###
#############################
            
            section = section + 1
            sec2text = sprintf("<hr><h2><a name = \"S2\">Section %s: Homogeneity between arrays</a></h2>", section)
            
            figure = figure + 1
            
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

            if(all(intgroup %in% names(phenoData(expressionset)@data)) && grouprep)
              {
                covar = pData(expressionset)[colnames(pData(expressionset))==intgroup[1]][,1]
                lev = levels(as.factor(covar))
                colourCovd = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[1])

                mplot2 = makePlot(con=con, name = "boxplot",
                  w=widthbox, h=8, fun = function() {
                    par(cex.axis = 1, pty = "s", lheight =((1/log10(numArrays))*long), mai = mai , omi = c(0,0,0,0), xaxt = xaxt)
                
                    if(expressionset@arrayInfo$channels == "single")
                      boxplotBeads(expressionset, las = 3, outline = FALSE,  col = colourCovd[as.factor(covar)], names = xname, ylab = "log(intensity)")
                
                    if(expressionset@arrayInfo$channels == "two")
                      {
                        nf = layout(matrix(1:3,1,3,byrow=TRUE),
                          c(2,2,2), 2, TRUE)
                        boxplotBeads(expressionset, las = 3, outline = FALSE, col = colourCovd[as.factor(covar)], names = xname, ylab = "log(red intensity)", whatToPlot = "R")
                        boxplotBeads(expressionset, las = 3, outline = FALSE, col = colourCovd[as.factor(covar)], names = xname, ylab = "log(green intensity)", whatToPlot = "G")              
                        boxplotBeads(expressionset, las = 3, outline = FALSE, col = colourCovd[as.factor(covar)], names = xname, ylab = "log(ratio)", whatToPlot = "M")
                      }}, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</center></b></td></table>\n", title = "Boxplots", fig = figure)
              } else {
                mplot2 = makePlot(con=con, name = "boxplot",
                  w=widthbox, h=8, fun = function() {
                    par(cex.axis = 1, pty = "s", lheight =((1/log10(numArrays))*long), mai = mai , omi = c(0,0,0,0), xaxt = xaxt)
                
                    if(expressionset@arrayInfo$channels == "single")
                      boxplotBeads(expressionset, las = 3, outline = FALSE, col= colours[2], names = xname, ylab = "log(intensity)")
                
                    if(expressionset@arrayInfo$channels == "two")
                      {
                        nf = layout(matrix(1:3,1,3,byrow=TRUE),
                          c(2,2,2), 2, TRUE)
                        boxplotBeads(expressionset, las = 3, outline = FALSE, col= colours[6], names = xname, ylab = "log(red intensity)", whatToPlot = "R")
                        boxplotBeads(expressionset, las = 3, outline = FALSE, col= colours[4], names = xname, ylab = "log(green intensity)", whatToPlot = "G")              
                        boxplotBeads(expressionset, las = 3, outline = FALSE, col= colours[2], names = xname, ylab = "log(ratio)", whatToPlot = "M")
                      }}, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.1\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><br><b>Figure %s</center></b></td></table>\n", title = "Boxplots", fig = figure)


              }
            htmltext2 = mplot2[[2]]

            if(expressionset@arrayInfo$channels == "single")
              boxspe = ""
            if(expressionset@arrayInfo$channels == "two")
              boxspe = " The left panel corresponds to the red channel. The middle panel shows the green channel. The right panel shows the boxplots of log<sub>2</sub>(ratio)."
           

            legendhom1 = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> presents boxplots of the log<sub>2</sub>(Intensities). Each box corresponds to one array. %s Typically, one expects the boxes to have similar size (IQR) and y position (median).</DIV>", figure, boxspe)          
      
############################
###Section 2.2 : Density ###
############################
            
            figure = figure + 1
            repe = ceiling(numArrays/9)
            colour = rep(brewer.pal(9, "Set1"),repe)

            if(all(intgroup %in% names(phenoData(expressionset)@data)) && grouprep)
              {
                covar = pData(expressionset)[colnames(pData(expressionset))==intgroup[1]][,1]
                lev = levels(as.factor(covar))
                colourCovd = brewer.pal(8,rownames(brewer.pal.info[brewer.pal.info$category=="qual",])[1])
                
                mplot3 = makePlot(con=con, name = "density",
                  w=10, h=10, fun = function() {
                    for(n in seq_len(max(group)))
                      {
                        if(expressionset@arrayInfo$channels == "single")
                          {
                            nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(2.8,2.8),c(1.8,2), TRUE)
                            par(xaxt = "n", mar = c(0,5,2,5))
                            plotBeadDensities(expressionset, what = "G", arrays = c(which(group==n)), col=colourCovd[as.factor(covar[group==n])], xlab="", ylab="", main="", xlim=xlim)
                            legend("topright",legend=lev,lty=1,col=colourCovd[as.factor(lev)], bty = "n", cex = 0.9)
                            mtext("Density", side = 2, adj = 0.5, padj = -4 , cex = 0.9)
                            if(max(group) != 1)
                              mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -1)
                            par(xaxt = "s",  mar = c(4,5,0,5))
                            multi("ecdf",d[group==n],xlim=xlim,"ECDF","log(intensity)","", col=colourCovd[as.factor(covar[group==n])])
                          }
                    
                        if(expressionset@arrayInfo$channels == "two")
                          {
                            nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),
                                         c(2.3,2,2), c(2,2), TRUE)
                            par(mar=c(0,5,1,1),xaxt = "n")
                            plotBeadDensities(expressionset, what = "R", arrays = c(which(group==n)), col =colourCovd[as.factor(covar[group==n])], xlab="", ylab="", main="")
                            legend("topright",legend=lev,lty=1,col=colourCovd[as.factor(lev)], bty = "n", cex = 0.9)
                            mtext("Density", side = 2, adj = 0.5, padj = -4 , cex = 0.9)
                            par(mar=c(0,2,1,1),xaxt = "n")
                            plotBeadDensities(expressionset, what = "G", arrays = c(which(group==n)), col =colourCovd[as.factor(covar[group==n])], xlab="", ylab="", main="")                   
                            legend("topright",legend=lev,lty=1,col=colourCovd[as.factor(lev)], bty = "n", cex = 0.9)
                            if(max(group) != 1)
                              mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -3 ,cex = 1)
                            par(mar=c(0,2,1,1),xaxt = "n")
                            plotBeadDensities(expressionset, what = "M", arrays = c(which(group==n)), col =colourCovd[as.factor(covar[group==n])], xlab="", ylab="", main="")
                            legend("topright",legend=lev,lty=1,col=colourCovd[as.factor(lev)], bty = "n", cex = 0.9)
                            par(mar=c(1,5,0,1), xaxt = "s")
                            multi("ecdf",dr[group==n],xlim=xlimr,"ECDF","log(red intensity)","", col =colourCovd[as.factor(covar[group==n])])
                            par(mar=c(1,2,0,1), xaxt = "s")
                            multi("ecdf",dg[group==n],xlim=xlimg,"","log(green intensity)","", col =colourCovd[as.factor(covar[group==n])])
                            par(mar=c(1,2,0,1), xaxt = "s")
                            multi("ecdf",dlr[group==n],xlim=xlim,"","log(ratio)","", col =colourCovd[as.factor(covar[group==n])])
                          }
                      }
                  }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", title="Density plots", fig = figure)
              } else {          
            
                mplot3 = makePlot(con=con, name = "density",
                  w=10, h=10, fun = function() {
                    for(n in seq_len(max(group)))
                      {
                        if(expressionset@arrayInfo$channels == "single")
                          {
                            nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(2.8,2.8),c(1.8,2), TRUE)
                            par(xaxt = "n", mar = c(0,5,2,5))
                            plotBeadDensities(expressionset, what = "G", arrays = c(which(group==n)), col = colour[seq_len(numArrays)], xlab="", ylab="", main="", xlim=xlim)
                            legend("topright", legend=sN[group==n],lty=1,col=brewer.pal(9, "Set1"), bty = "n", cex = 0.9)
                            mtext("Density", side = 2, adj = 0.5, padj = -4 , cex = 0.9)
                            if(max(group) != 1)
                              mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -1)
                            par(xaxt = "s",  mar = c(4,5,0,5))
                            multi("ecdf",d[group==n],xlim=xlim,"ECDF","log(intensity)","")
                          }
                    
                        if(expressionset@arrayInfo$channels == "two")
                          {
                            nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),
                                         c(2.3,2,2), c(2,2), TRUE)
                            par(mar=c(0,5,1,1),xaxt = "n")
                            plotBeadDensities(expressionset, what = "R", arrays = c(which(group==n)), col = colour[seq_len(numArrays)], xlab="", ylab="", main="")
                            legend("topright",legend=sN[group==n],lty=1,col=brewer.pal(9, "Set1"), bty = "n", cex = 0.9)
                            mtext("Density", side = 2, adj = 0.5, padj = -4 , cex = 0.9)
                            par(mar=c(0,2,1,1),xaxt = "n")
                            plotBeadDensities(expressionset, what = "G", arrays = c(which(group==n)), col = colour[seq_len(numArrays)], xlab="", ylab="", main="")                   
                            legend("topright",legend= sN[group==n],lty=1,col=brewer.pal(9, "Set1"), bty = "n", cex = 0.9)
                            if(max(group) != 1)
                              mtext(paste("Group ", n, sep = ""),side = 3,adj = 0.5, padj = -3 ,cex = 1)
                            par(mar=c(0,2,1,1),xaxt = "n")
                            plotBeadDensities(expressionset, what = "M", arrays = c(which(group==n)), col = colour[seq_len(numArrays)], xlab="", ylab="", main="")                   
                            legend("topright", legend=sN[group==n],lty=1,col=brewer.pal(9, "Set1"), bty = "n", cex = 0.9)
                            par(mar=c(1,5,0,1), xaxt = "s")
                            multi("ecdf",dr[group==n],xlim=xlimr,"ECDF","log(red intensity)","")
                            par(mar=c(1,2,0,1), xaxt = "s")
                            multi("ecdf",dg[group==n],xlim=xlimg,"","log(green intensity)","")
                            par(mar=c(1,2,0,1), xaxt = "s")
                            multi("ecdf",dlr[group==n],xlim=xlim,"","log(ratio)","")
                          }
                      }
                  }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><center><a name = \"S2.2\"><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></a></A><center><br><b>Figure %s</b></center></td></table>\n", title="Density plots", fig = figure)
              }
            htmltext3 = mplot3[[2]]            
     
            legendhom2 = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows density estimates (histograms) of the data. Arrays whose distributions are very different from the others should be considered for possible problems.</DIV>", figure)          
           
##########################################
###Section 4 : Between array comparison###
##########################################
            
            ##Heatmap
            hmap = hmap(expressionset, sN, section, figure, outM, numArrays, intgroup)
            
            section = hmap$section
            figure = hmap$figure
            htmltext4 = hmap$htmltext4
            sec4text = hmap$sec4text
            legendheatmap = hmap$legendheatmap
            
##########################################
###Section 5 : Variance Mean Dependence###
##########################################
            
            ##meanSdplot
            msdp = msdp(expressionset, section, figure, con, dat)
            
            section = msdp$section
            figure = msdp$figure
            htmltext5 = msdp$htmltext5
            sec5text = msdp$sec5text
            legendsdmean = msdp$legendsdmean

##########################
### Scores computation ###
##########################
            sc = scores(expressionset=expressionset,
              numArrays=numArrays, M=M, ldat=NULL, outM=outM, dat=dat,
              rc=NULL, gc=NULL, maxc=NULL, maxr=NULL, nuse=NULL,
              rle=NULL) 
            
            scores = matrix("",ncol=length(sc),nrow=numArrays)
            for( i in seq_len(length(sc)))
              scores[unlist(sc[[i]][seq_len(length(sc[[i]]))]),i]="*"



##########################
### Writing the report ###
##########################

            arg = as.list(match.call(expand.dots = TRUE))
            
            con = report(expressionset=expressionset, arg=arg,
            sNt=sNt, sN=sN, sec1text=sec1text, mapdf=mapdf,
            matext1=matext1, nfig=nfig, legendMA=legendMA,
            batext=batext, nfig2=nfig2, bapng=bapng, ftext=ftext,
            pttext=NULL, legendpt=NULL, nfig3=nfig3, fpng=fpng,
            legendlocal=legendlocal, sec2text=sec2text,
            htmltext2=htmltext2, legendhom1=legendhom1, group=group,
            htmltext3=htmltext3, legendhom2=legendhom2,
            sec4text=sec4text,
            htmltext4=htmltext4, legendheatmap=legendheatmap,
            sec5text=sec5text, htmltext5=htmltext5,
            legendsdmean=legendsdmean, scores=scores) 
            writeLines("</table>", con)

            z=sessionInfo("arrayQualityMetrics")
            version = z$otherPkgs[[1]]$Version
            rversion = sessionInfo()$R.version$version.string
            writeLines(sprintf("<hr><DIV style=\"font-size: 13; font-family: Lucida Grande\">This report has been created with arrayQualityMetrics %s under %s</DIV>",version, rversion), con)

            closeHtmlPage(con)
            
          } ####end set method BeadLevelList
          )
