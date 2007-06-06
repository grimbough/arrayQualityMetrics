setGeneric("arrayQualityMetrics", function(expressionset, outfile, numberofgraphs) standardGeneric("arrayQualityMetrics"))

##library("Biobase")
##library("affyPLM")
##library("RColorBrewer")
##library("geneplotter")
##library("arrayMagic")

##to be replaced by loading genefilter when it works
dist2 = function(x, fun=function(a,b) mad(a-b)) {
  res = matrix(as.numeric(NA), ncol=ncol(x), nrow=ncol(x))
  colnames(res) = rownames(res) = colnames(x)
  if(ncol(x)>=2) {
    for(j in 2:ncol(x))
      for(i in 1:(j-1))
        res[i, j] = res[j, i] = fun(x[,i], x[,j])
  }
  return(res)
}

##lists
mat2list = function(x)
      lapply(seq_len(ncol(x)), function(i) x[,i])

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


##foreground and background representations
localisation = function(x, title1, title2)
  {
    plot.imageMatrix(x,
                     xlab = "",
                     ylab = "",
                     zScale = FALSE,
                     width = 4,
                     height = 7)
    mtext(title1,side = 3,adj = 0.5, padj = -1 ,cex = 0.7)
    mtext(title2,side = 2,adj = 0.5, padj = -1 ,cex = 0.7)
  }

##makePlot function
makePlot = function(con, name, w, h=devDims(w)$height, fun, psz=12, print=TRUE,  isPlatePlot=FALSE, isImageScreen=FALSE) {

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

  if (print)
    cat(sprintf("<CENTER><A HREF=\"%s\"><IMG SRC=\"%s\"/></A></CENTER><BR>\n",
                outf[1], outf[2]), file=con)

  return(res)
}

##NChannelSet

setMethod("arrayQualityMetrics",signature(expressionset ="NChannelSet"),
function(expressionset, outfile, numberofgraphs)
  {
      fn  = paste(outfile,"_QMreport",sep="")
      title = paste(outfile, " quality metrics report", sep="")
      titletext = sprintf("<hr><h1><center>%s</h1></center><hr><table border = \"0\" cellspacing = 5 cellpadding = 2><CENTER>", title)
      con = openHtmlPage(fn, title)
      writeLines(titletext, con)
      
      if(all(assayData(expressionset)$R >= 0) == FALSE)
        rc = assayData(expressionset)$R
      if(all(assayData(expressionset)$R >= 0) == TRUE)
        rc = log2(assayData(expressionset)$R)
      
      if(all(assayData(expressionset)$G >= 0) == FALSE)
        gc = assayData(expressionset)$G
      if(all(assayData(expressionset)$G >= 0) == TRUE)
        gc = log2(assayData(expressionset)$G)
   
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
      if(length(phenoData(expressionset)$dyeswap) != 0)
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
      if(exists("sNt"))
        {
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
      writeLines("</table><hr>", con)
    }
      
      ##Heatmap
      colourRange = brewer.pal(9,"Greys")
      outM = as.dist(dist2(na.omit(dat)))
      marg = ceiling(long/(sqrt(long)/1.5)+log10(numArrays))
      
      makePlot(con=con, name = paste(outfile, "_heatmap", sep = ""),
                w=4, h=4, fun = function() {
       heatmap(as.matrix(outM),
              labRow = sN,
              labCol = sN,
              col = colourRange,
              scale = "none",
              main = "",
              margins = c(marg,marg))
              }, print=TRUE)

      
      
      hpng = paste(outfile, "_heatmap.png", sep = "")  
      hpdf = paste(outfile, "_heatmap.pdf", sep = "")
      png(file= hpng)
      heatmap(as.matrix(outM),
              labRow = sN,
              labCol = sN,
              col = colourRange,
              scale = "none",
              main = "",
              margins = c(marg,marg))
      dev.copy(pdf, file = hpdf)
      dev.off()
      dev.off()
      
      m = matrix(pretty(outM,9),nrow=1,ncol=length(pretty(outM,9)))
      hlpng = paste(outfile, "_heatmaplegend.png", sep = "")  
      png(file= hlpng, width = 150, height = 450)
      image(m,xaxt="n",yaxt="n",ylab="Distance", col = colourRange, cex.lab = 0.8, mgp = c(1.5,1,0) )
      axis(2, label= as.list(pretty(outM,9)),at=seq(0,1,by=(1/(length(pretty(outM,9))-1))), cex.axis = 0.8, padj = 1)
      dev.off()
      htext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></CENTER><BR></td></tr></table><hr>\n", "Heatmap representation of the distance between experiments", hpdf, hpng, hlpng)
      writeLines(htext, con) 
      
      hc = hclust(outM)
      
      if(!missing(numberofgraphs))
        k = numberofgraphs
      if(missing(numberofgraphs) && numArrays > 50)
        k = ceiling(numArrays/50)
      if(missing(numberofgraphs) && numArrays <= 50)
        k = 1    
      group = cutree(hc, k = k)
      
      
      ##Boxplots
      colours = brewer.pal(12, "Paired")
      bpng = paste(outfile, "_boxplot.png", sep = "")
      bpdf = paste(outfile, "_boxplot.pdf", sep = "")
      png(file = bpng)
      col = c(rep(colours[4],numArrays),rep(colours[6],numArrays))
      if(numArrays <= 15)
        {
          xname1 = c(sN,sN)
          xname2 = sN
          mai = c(0,0.4,0.2,0.2)
          xaxt = "s"
        }
      if(numArrays > 15 && numArrays <= 50)
        {
          xname1 = FALSE
          xname2 = sN
          mai = c((log(long)-log(long)/1.2),0.4,0.2,0.2)
          xaxt = "s"
        }
      if(numArrays > 50)
        {
          xname1 = FALSE
          xname2 = FALSE
          mai = c(0,0.4,0.2,0.2)
          xaxt = "n"
        }
      
      par(mfrow = c(1,2),
          cex.axis = (1/log10(numArrays)-0.1),
          pty = "s",
          lheight =((1/log10(numArrays))*long),
          mai = mai,
          omi = c(0,0,0,0),
          xaxt = xaxt)
      boxplot(c(lgreenc, lredc), col = col, las = 3, range = 0,
              names = xname1)
      boxplot(ldat, col = colours[2], las = 3, range = 0,
              names = xname2)
      dev.copy(pdf, file = bpdf)
      dev.off()
      dev.off()
      btext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></td></tr></table><hr>\n", "Boxplots", bpdf, bpng)
      writeLines(btext, con) 
      
      ##Density if 1 group
      if(max(group) == 1)
        {
          dpng = paste(outfile,"_density.png", sep = "")
          dpdf = paste(outfile,"_density.pdf", sep = "")
          png(file = dpng)
          xlim = c(min(na.omit(dat)),max(na.omit(dat)))
          xlimr = c(min(na.omit(rc)),max(na.omit(rc)))
          xlimg = c(min(na.omit(gc)),max(na.omit(gc)))
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
          multi("ecdf",ldat,xlim,"","log(ratio)","")
          dev.copy(pdf, file = dpdf)
          dev.off()
          dev.off()
        }
      
      ##Density if more than 1 group
      if(max(group) > 1)
        {
          dpng = paste(outfile,"_density.png", sep = "")
          dpdf = paste(outfile,"_density.pdf", sep = "")
          pdf(file = dpdf)
          xlim = c(min(na.omit(dat)),max(na.omit(dat)))
          xlimr = c(min(na.omit(rc)),max(na.omit(rc)))
          xlimg = c(min(na.omit(gc)),max(na.omit(gc)))
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
        }
      dtext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "Density plots", dpdf, dpng)
      writeLines(dtext, con)
      
      ##Background rank representation
      
      if("r" %in% rownames(featureData(expressionset)@varMetadata) && "c" %in% rownames(featureData(expressionset)@varMetadata))
        {
          r = featureData(expressionset)$r
          c = featureData(expressionset)$c
          
          if("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))
            {
              bapng = paste(outfile, "_background.png", sep = "")
              bapdf = paste(outfile, "_background.pdf", sep = "")
              pdf(file = bapdf)
              nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = FALSE),
                           c(1.2,1.2,1.2), c(2,2), TRUE)
              for(a in 1:numArrays)
                {
                  intr = cbind(r,c,rbg[,a])
                  intg = cbind(r,c,gbg[,a])
                  
                  re = matrix(NA,ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
                  g = matrix(NA,ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
                  
                  for(i in 1:nrow(intr))
                    {
                      re[intr[i,1],intr[i,2]] = intr[i,3]
                      g[intg[i,1],intg[i,2]] = intg[i,3]
                    }
                  
                  mr = rank(re)
                  mg = rank(g)
                  par(xaxt = "n", yaxt = "n",mar=c(2,2,2,1))
                  localisation(mr,sN[a],"Red Background")
                  localisation(mg,sN[a],"Green Background")
                }   
              dev.off()
              png(file = bapng, width = 200, height = 200)
              nf <- layout(matrix(c(1,2),1,2,byrow = FALSE),
                           c(1.2,1.2), c(2,2), TRUE)
              par(xaxt = "n", yaxt = "n",mar=c(2,2,2,1))
              localisation(mr,"","Red Intensity")
              localisation(mg,"","Green Intensity")
              dev.off()
            }
          batext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "Background representation on the array", bapdf, bapng)
          writeLines(batext, con)
          
          ##Foreground rank representation
          
          fpng = paste(outfile, "_foreground.png", sep = "")
          fpdf = paste(outfile, "_foreground.pdf", sep = "")
          pdf(file = fpdf)
          nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = FALSE),
                       c(1.2,1.2,1.2), c(2,2), TRUE)
          for(a in 1:numArrays)
            {
              intrf = cbind(r,c,rc[,a])
              intgf = cbind(r,c,gc[,a])
              
              rf = matrix(NA,ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
              gf = matrix(NA,ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
              
              for(i in 1:nrow(intrf))
                {
                  rf[intrf[i,1],intrf[i,2]] = intrf[i,3]
                  gf[intgf[i,1],intgf[i,2]] = intgf[i,3]
                }
                
              mrf = rank(rf)
              mgf = rank(gf)
              par(xaxt = "n", yaxt = "n",mar=c(2,2,2,1))
              localisation(mrf,sN[a],"Red Intensity")
              localisation(mgf,sN[a],"Green Intensity")
            }
          dev.off()
          png(file = fpng, width = 200, height = 200)
          nf <- layout(matrix(c(1,2),1,2,byrow = FALSE),
                       c(1.2,1.2), c(2,2), TRUE)
          par(xaxt = "n", yaxt = "n",mar=c(2,2,2,1))
          localisation(mrf,"","Red Intensity")
          localisation(mgf,"","Green Intensity")
          dev.off()
          ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "Foreground representation on the array", fpdf, fpng)
          writeLines(ftext, con)
        }

      ##GC
    
      if("GC" %in% rownames(featureData(expressionset)@varMetadata))
        {
          gcpdf = paste(outfile,"_GCcontent.pdf", sep = "")
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
        
          gcopng = paste(outfile,"_overall_GCcontent.png", sep = "")
          gcopdf = paste(outfile,"_overall_GCcontent.pdf", sep = "")
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
          
          gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "GC content effect ", gcpdf, "per array", " and ", gcopdf, "global", gcopdf, gcopng)
          writeLines(gctext, con)
        }
      
      ##Scatter plots between replicates
      ##Have to find a better solution for the xlim and ylim.
      if("replicates" %in% names(phenoData(expressionset)))
        {
          rpng = paste(outfile,"_replicates.png", sep = "")
          rpdf = paste(outfile,"_replicates.pdf", sep = "")
          pdf(file = rpdf)
          lev = levels(expressionset@phenoData$replicates)
          for(i in 1:length(lev))
            {
              subset = dat[, names(phenoData(expressionset)$replicates[
                  phenoData(expressionset)$replicates == lev[i]])]
              comp = combn(ncol(subset),2)
              layout(matrix(1:16,4,4,byrow = TRUE), rep(2,4), rep(2,4), TRUE)
              par(mar = c(0,0,0,0))
              for(j in 1:ncol(comp))
                {
                  smoothScatter(subset[,comp[1,j]],subset[,comp[2,j]], xaxt = "n", xlim = xlim, yaxt = "n", ylim = c(min(dat),max(dat)), xlab = "", ylab = "")
                  abline(a=0,b=1,col="red")
                }
            }
          dev.copy(png, file = rpng)
          dev.off()
          dev.off()
          rtext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "Scatter plots between replicates", rpdf, rpng)
          writeLines(rtext, con)
        }
      
      ##MAplots
      ##function from affyQCReport
      M =  rc - gc
      A = 0.5*(rc + gc)
      
      app = 4 + 2*(sum(numArrays>c(4,6)))
      nfig = ceiling(numArrays/8)
      
      plotNames = paste(outfile,"_MA", 1:nfig, sep="")
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
    
      matext1 = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></td><td><A HREF=\"%s\">%s</td></A></CENTER><BR></tr></td>\n", "MvA plots",  mapdf[1],  mapng[1],mapdf[1], "MvA plot 1")
      writeLines(matext1, con)
      if(nfig >= 2)
        {
          for(i in 2:nfig)
            {
              matext2 = sprintf("<tr><td></td><td></td><td><A HREF=\"%s\">%s%s</A></CENTER><BR></tr></td>\n", mapdf[i], "MvA plot ", i)
          
              writeLines(matext2, con)
            }
        }
      writeLines("</table><hr>", con)
      

      ##Mapping of probes
    
      if("gn" %in% rownames(featureData(expressionset)@varMetadata))
        {
          probenames = expressionset@featureData$gn
          genes = grep("^NM_", probenames)
          facgene = as.vector(probenames)
          
          facgene[genes] = 1
          facgene[-genes] = 0
          gpdf = paste(outfile,"_GenesMapping.pdf", sep = "")
          gopng = paste(outfile,"_overall_GenesMapping.png", sep = "")
          gopdf = paste(outfile,"_overall_GenesMapping.pdf", sep = "")
          
          cols = brewer.pal(9, "Set1")
          pdf(gpdf)
          nf <- layout(matrix(1:16,4,4,byrow=TRUE), c(2,1.8,1.8,1.8), c(1.8,1.8,1.8,2), FALSE)
          for(a in 1:numArrays)
            {
              if(a %in% c(seq(14,numArrays,by=16),seq(15,numArrays,by=16),seq(16,numArrays,by=16)))
                {
                  par(mar = c(2,0,0,0))
                  multi("density",dat[,a]~facgene,xlim,"","","", col = cols[c(9,2)], yaxt = "n", ylim = c(0,1))
                  legend("topright",legend=sN[a], bty = "n", cex = 0.6)
                }
              if(a %in%  c(seq(1,numArrays,by=16),seq(5,numArrays,by=16),seq(9,numArrays,by=16)))
                {
                  par(mar = c(0,2,0,0))
                  multi("density",dat[,a]~facgene,xlim,"","","", col = cols[c(9,2)], xaxt = "n", ylim = c(0,1))
                  legend("topright",legend=sN[a], bty = "n", cex = 0.6)
                }
              if(a %in% seq(13,numArrays,by=16))
                {
                  par(mar = c(2,2,0,0))
                  multi("density",dat[,a]~facgene,xlim,"","","", col = cols[c(9,2)], ylim = c(0,1))
                  legend("topright",legend=sN[a], bty = "n", cex = 0.6)
                }
              if(!(a %in%  seq(13,numArrays,by=16)) && !(a %in%  c(seq(1,numArrays,by=16),seq(5,numArrays,by=16),seq(9,numArrays,by=16))) && !(a %in% c(seq(14,numArrays,by=16),seq(15,numArrays,by=16),seq(16,numArrays,by=16))))
                {
                  par(mar = c(0,0,0,0))
                  multi("density",dat[,a]~facgene,xlim,"","","", col = cols[c(9,2)],  xaxt = "n", yaxt = "n",ylim = c(0,1))
                  legend("topright",legend=sN[a], bty = "n", cex = 0.6)
                }
            }
          dev.off()
    
          png(file = gopng)
          multi("density",dat~facgene,xlim,"","","", col = cols[c(9,2)])
          legend("topright", c("Unigene","Unmapped"),lty=1,lwd=2,col= c(cols[c(2,9)]), bty = "n")
          dev.copy(pdf, file = gopdf)
          dev.off()
          dev.off()

          gotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "Gene mapping ", gpdf, "per array", " and ", gopdf, "global", gopdf, gopng)
          writeLines(gotext, con)

        }
      writeLines("</table>", con)
      closeHtmlPage(con)
    
    }#endofarrayQualityMetrics
)

##plotting functions for ExpressionSet

aqm.expressionset = function(expressionset, outfile, numberofgraphs, con)
  {
    ##data preparation
    dat = if(all(exprs(expressionset) < 0, na.rm=TRUE)) {
       exprs(expressionset)
     } else {
       log2(exprs(expressionset))
     }

    
    sN = sampleNames(expressionset)
    ##list to use multidensity, multiecdf and boxplot
    ldat = mat2list(dat)
    gN = featureNames(expressionset)
    numArrays = ncol(dat)
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
    if(exists("sNt"))
      {
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
        writeLines("</table><hr>", con)
      }
                
    ##Heatmap
    colourRange = brewer.pal(9,"Greys")
    outM = as.dist(dist2(na.omit(dat)))
    marg = ceiling(long/(sqrt(long)/1.5)+log10(numArrays))
                
    hpng = paste(outfile, "_heatmap.png", sep = "")  
    hpdf = paste(outfile, "_heatmap.pdf", sep = "")
    png(file= hpng)
    heatmap(as.matrix(outM),
            labRow = sN,
            labCol = sN,
            col = colourRange,
            scale = "none",
            main = "",
            margins = c(marg,marg))
    dev.copy(pdf, file = hpdf)
    dev.off()
    dev.off()
    
    m = matrix(pretty(outM,9),nrow=1,ncol=length(pretty(outM,9)))
    hlpng = paste(outfile, "_heatmaplegend.png", sep = "")  
    png(file= hlpng, width = 150, height = 450)
    image(m,xaxt="n",yaxt="n",ylab="Distance", col = colourRange, cex.lab = 0.8, mgp = c(1.5,1,0) )
    axis(2, label= as.list(pretty(outM,9)),at=seq(0,1,by=(1/(length(pretty(outM,9))-1))), cex.axis = 0.8, padj = 1)
    dev.off()
    htext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></CENTER><BR></td></tr></table><hr>\n", "Heatmap representation of the distance between experiments", hpdf, hpng, hlpng)
    writeLines(htext, con) 
    
    hc = hclust(outM)
    
    if(!missing(numberofgraphs))
      k = numberofgraphs
    if(missing(numberofgraphs) && numArrays > 50)
      k = ceiling(numArrays/50)
    if(missing(numberofgraphs) && numArrays <= 50)
      k = 1    
    group = cutree(hc, k = k)
                
    ##Boxplots
    colours = brewer.pal(12, "Paired")
    bpng = paste(outfile, "_boxplot.png", sep = "")
    bpdf = paste(outfile, "_boxplot.pdf", sep = "")
    png(file = bpng)
    par(cex.axis = 0.8, pty = "s", lheight =((1/log10(numArrays))*long), mai = c(((long/12)+0.2),0.4,0.2,0.2) , omi = c(0,0,0,0))
    boxplot(ldat, col = colours[2], las = 3, range = 0, names = sN)
    dev.copy(pdf, file = bpdf)
    dev.off()
    dev.off()
    btext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></td></tr></table><hr>\n", "Boxplots", bpdf, bpng)
    writeLines(btext, con) 
    
    ##Density if 1 group
    dpng = paste(outfile,"_density.png", sep = "")
    dpdf = paste(outfile,"_density.pdf", sep = "")
    png(file = dpng)
    xlim = c(min(na.omit(dat)),max(na.omit(dat)))
    nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(2.8,2.8),c(1.8,2), TRUE)
    par(xaxt = "n", cex.axis = 0.8, mar = c(0,5,2,5))
    multi("density",ldat,xlim,"Density","","")
    par(xaxt = "s", cex.axis = 0.8, mar = c(4,5,0,5))
    multi("ecdf",ldat,xlim,"ECDF","log(intensity)","")
    dev.copy(pdf, file = dpdf)
    dev.off()
    dev.off()
                
    ##Density if more than 1 group
    dpng = paste(outfile,"_density.png", sep = "")
    dpdf = paste(outfile,"_density.pdf", sep = "")
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
    dtext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "Density plots", dpdf, dpng)
    writeLines(dtext, con)
    
    ##GC
    if("GC" %in% rownames(featureData(expressionset)@varMetadata))
      {
        gcpdf = paste(outfile,"_GCcontent.pdf", sep = "")
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
        
        gcopng = paste(outfile,"_overall_GCcontent.png", sep = "")
        gcopdf = paste(outfile,"_overall_GCcontent.pdf", sep = "")
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
                
        gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "GC content effect ", gcpdf, "per array", " and ", gcopdf, "global", gcopdf, gcopng)
        writeLines(gctext, con)
      }
                
    ##Scatter plots between replicates
    ##Have to find a better solution for the xlim and ylim.
    if("replicates" %in% names(phenoData(expressionset)))
      {
        rpng = paste(outfile,"_replicates.png", sep = "")
        rpdf = paste(outfile,"_replicates.pdf", sep = "")
        pdf(file = rpdf)
        lev = levels(expressionset@phenoData$replicates)
        for(i in 1:length(lev))
          {
            subset = dat[, names(phenoData(expressionset)$replicates[
                phenoData(expressionset)$replicates == lev[i]])]
            comp = combn(ncol(subset),2)
            layout(matrix(1:16,4,4,byrow = TRUE), rep(2,4), rep(2,4), TRUE)
            par(mar = c(0,0,0,0))
            for(j in 1:ncol(comp))
              {
                smoothScatter(subset[,comp[1,j]],subset[,comp[2,j]], xaxt = "n", xlim = xlim, yaxt = "n", ylim = c(min(dat),max(dat)), xlab = "", ylab = "")
                abline(a=0,b=1,col="red")
              }
          }
        dev.copy(png, file = rpng)
        dev.off()
        dev.off()
        rtext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "Scatter plots between replicates", rpdf, rpng)
        writeLines(rtext, con)
      }
    
    ##MAplots
    ##function from affyQCReport
    medArray = rowMedians(dat)
    M =  dat-medArray
    A = (dat+medArray)/2
    app = 4 + 2*(sum(numArrays>c(4,6)))
    nfig = ceiling(numArrays/8)
    
    plotNames = paste(outfile,"_MA", 1:nfig, sep="")
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
    
    matext1 = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></td><td><A HREF=\"%s\">%s</td></A></CENTER><BR></tr></td>\n", "MvA plots",  mapdf[1],  mapng[1],mapdf[1], "MvA plot 1")
    writeLines(matext1, con)
    if(nfig >= 2)
      {
        for(i in 2:nfig)
          {
            matext2 = sprintf("<tr><td></td><td></td><td><A HREF=\"%s\">%s%s</A></CENTER><BR></tr></td>\n", mapdf[i], "MvA plot ", i)
            writeLines(matext2, con)
          }
      }
    writeLines("</table><hr>", con)
                
    ##Mapping of probes
    if("gn" %in% rownames(featureData(expressionset)@varMetadata))
      {
        probenames = expressionset@featureData$gn
        genes = grep("^NM_", probenames)
        facgene = as.vector(probenames)
                    
        facgene[genes] = 1
        facgene[-genes] = 0
        gpdf = paste(outfile,"_GenesMapping.pdf", sep = "")
        gopng = paste(outfile,"_overall_GenesMapping.png", sep = "")
        gopdf = paste(outfile,"_overall_GenesMapping.pdf", sep = "")
        
        cols = brewer.pal(9, "Set1")
        pdf(gpdf)
        nf <- layout(matrix(1:16,4,4,byrow=TRUE), c(2,1.8,1.8,1.8), c(1.8,1.8,1.8,2), FALSE)
        for(a in 1:numArrays)
          {
            if(a %in% c(seq(14,numArrays,by=16),seq(15,numArrays,by=16),seq(16,numArrays,by=16)))
              {
                par(mar = c(2,0,0,0))
                multi("density",dat[,a]~facgene,xlim,"","","", col = cols[c(9,2)], yaxt = "n", ylim = c(0,1))
                legend("topright",legend=sN[a], bty = "n", cex = 0.6)
              }
            if(a %in%  c(seq(1,numArrays,by=16),seq(5,numArrays,by=16),seq(9,numArrays,by=16)))
              {
                par(mar = c(0,2,0,0))
                multi("density",dat[,a]~facgene,xlim,"","","", col = cols[c(9,2)], xaxt = "n", ylim = c(0,1))
                legend("topright",legend=sN[a], bty = "n", cex = 0.6)
              }
            if(a %in% seq(13,numArrays,by=16))
              {
                par(mar = c(2,2,0,0))
                multi("density",dat[,a]~facgene,xlim,"","","", col = cols[c(9,2)], ylim = c(0,1))
                legend("topright",legend=sN[a], bty = "n", cex = 0.6)
              }
            if(!(a %in%  seq(13,numArrays,by=16)) && !(a %in%  c(seq(1,numArrays,by=16),seq(5,numArrays,by=16),seq(9,numArrays,by=16))) && !(a %in% c(seq(14,numArrays,by=16),seq(15,numArrays,by=16),seq(16,numArrays,by=16))))
              {
                par(mar = c(0,0,0,0))
                multi("density",dat[,a]~facgene,xlim,"","","", col = cols[c(9,2)],  xaxt = "n", yaxt = "n",ylim = c(0,1))
                legend("topright",legend=sN[a], bty = "n", cex = 0.6)
              }
          }
        dev.off()
        
        png(file = gopng)
        multi("density",dat~facgene,xlim,"","","", col = cols[c(9,2)])
        legend("topright", c("Unigene","Unmapped"),lty=1,lwd=2,col= c(cols[c(2,9)]), bty = "n")
        dev.copy(pdf, file = gopdf)
        dev.off()
        dev.off()
        
        gotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "Gene mapping ", gpdf, "per array", " and ", gopdf, "global", gopdf, gopng)
        writeLines(gotext, con)
      }
    l = list(numArrays,sN)
    return(l) 
  }

##ExpressionSet
setMethod("arrayQualityMetrics",signature(expressionset="ExpressionSet"),
          function(expressionset, outfile, numberofgraphs)
          {
            fn  = paste(outfile,"_QMreport",sep="")
            title = paste(outfile, " quality metrics report", sep="")
            titletext = sprintf("<hr><h1><center>%s</h1></center><hr><table border = \"0\" cellspacing = 5 cellpadding = 2><CENTER>", title)
            con = openHtmlPage(fn, title)
            writeLines(titletext, con)

            aqm.expressionset(expressionset, outfile, numberofgraphs, con)
               
            writeLines("</table>", con)
            closeHtmlPage(con)
          })#end set method ExpressionSet



##AffyBatch

setMethod("arrayQualityMetrics",signature(expressionset="AffyBatch"),
          function(expressionset, outfile, numberofgraphs){
            fn  = paste(outfile,"_QMreport",sep="")
            title = paste(outfile, " quality metrics report", sep="")
            titletext = sprintf("<hr><h1><center>%s</h1></center><hr><table border = \"0\" cellspacing = 5 cellpadding = 2><CENTER>", title)
            con = openHtmlPage(fn, title)
            writeLines(titletext, con)
            
            l = aqm.expressionset(expressionset, outfile, numberofgraphs, con)
            numArrays = as.numeric(l[1])
            sN = l[[2]]
              
            acol = sample(brewer.pal(8, "Dark2"), numArrays, replace = (8<numArrays))
            rnaDeg = AffyRNAdeg(expressionset)
            affypng1 = paste(outfile,"_RNAdeg.png", sep = "")
            affypdf1 = paste(outfile,"_RNAdeg.pdf", sep = "")
            png(file = affypng1)
            plotAffyRNAdeg(rnaDeg, cols = acol, lwd = 2)
            dev.copy(pdf, file = affypdf1)
            dev.off()
            dev.off()
                
            pp1 = preprocess(expressionset)
            dataPLM = fitPLM(pp1, background = FALSE, normalize = FALSE)
            affypng2 = paste(outfile,"_RLE.png", sep = "")
            affypdf2 = paste(outfile,"_RLE.pdf", sep = "")
            png(file = affypng2)
            Mbox(dataPLM, ylim = c(-1, 1), names = sN, col = "lightblue",
                 whisklty = 0, staplelty = 0, main = "RLE", las = 3)
            dev.copy(pdf, file = affypdf2)
            dev.off()
            dev.off()
            
            affypng3 = paste(outfile,"_NUSE.png", sep = "")
            affypdf3 = paste(outfile,"_NUSE.pdf", sep = "")
            png(file = affypng3)
            boxplot(dataPLM, ylim = c(0.95, 1.5), names = sN,
                    outline = FALSE, col = "lightblue", main = "NUSE", las = 2)
            dev.copy(pdf, file = affypdf3)
            dev.off()
            dev.off()
            affytext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER><BR></tr></td></table><hr>\n", "RNA degradation plot", affypdf1, affypng1, "RLE plot", affypdf2, affypng2, "NUSE plot", affypdf3, affypng3)
            writeLines(affytext, con)
            
            writeLines("</table>", con)
            closeHtmlPage(con)
          })#end set method AffyBatch

