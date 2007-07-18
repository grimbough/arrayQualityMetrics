setGeneric("arrayQualityMetrics", function(expressionset, outfile, log.transformed, numberofgraphs) standardGeneric("arrayQualityMetrics"))


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
  
  htmltext = sprintf(text, title, outf[1], outf[2], fig)
  
  writeLines(htmltext, con)
  
  return(res)
}

##NChannelSet

setMethod("arrayQualityMetrics",signature(expressionset = "NChannelSet"),
          function(expressionset, outfile, log.transformed, numberofgraphs)
          {
            fn  = paste(outfile,"_QMreport",sep="")
            title = paste(outfile, " quality metrics report", sep="")
            titletext = sprintf("<hr><h1><center>%s</h1></center><table border = \"0\" cellspacing = 5 cellpadding = 2>", title)
            con = openHtmlPage(fn, title)
            writeLines(titletext, con)
      
            ##data preparation
            if(log.transformed == FALSE)
              {
                rc = log2(assayData(expressionset)$R)
                gc = log2(assayData(expressionset)$G)
              }
            if(log.transformed == TRUE)
              {
                rc = assayData(expressionset)$R
                gc = assayData(expressionset)$G
              }
            
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

            outM = as.dist(dist2(na.omit(dat)))
            hc = hclust(outM)
      
            if(!missing(numberofgraphs))
              k = numberofgraphs
            if(missing(numberofgraphs) && numArrays > 50)
              k = ceiling(numArrays/50)
            if(missing(numberofgraphs) && numArrays <= 50)
              k = 1    
            group = cutree(hc, k = k)
     

##########################################
###Section 1 : Individual array quality###
##########################################
            
            section = 1
            figure = 1
            sec1text = sprintf("<hr><h2>Section %s: Individual array quality</h2>", section)
            writeLines(sec1text, con)
            
           ##MAplots
            ##function from affyQCReport
            M = rc - gc
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
    
            matext1 = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></td><td><A HREF=\"%s\">%s</td></A></tr></td>\n", "MvA plots",  mapdf[1],  mapng[1], figure , mapdf[1], "MvA plot 1")
            writeLines(matext1, con)
            if(nfig >= 2)
              {
                for(i in 2:nfig)
                  {
                    matext2 = sprintf("<tr><td></td><td></td><td><A HREF=\"%s\">%s%s</A><BR></tr></td>\n", mapdf[i], "MvA plot ", i)
          
                    writeLines(matext2, con)
                  }
              }
            writeLines("</table>", con)
            legendMA = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> represents MA plot for each array. <br> M A-plots are useful for pairwise comparisons between arrays. M and A are defined as :<br>
M = log2(I1) - log2(I2)<br>
A = 1/2 log2(I1)+log2(I2))<br>
where I1 and I2 are the vectors of normalized intensities of two channels, on the original data scale (i. e. not logarithm-transformed).Typically, we expect the mass of the distribution in an M A-plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A.
Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).</DIV>", figure)          
            writeLines(legendMA, con)
            
#################################
###Section 2 : Spatial Effects###
#################################

            ##Background rank representation
      
            if("r" %in% rownames(featureData(expressionset)@varMetadata) && "c" %in% rownames(featureData(expressionset)@varMetadata))
              {
                section = section + 1
                sec2text = sprintf("<hr><h2>Section %s: Spatial Effects</h2>", section)
                writeLines(sec2text, con)

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
                batext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER></tr></td></table>\n", "Background representation on the array", bapdf, bapng)
                writeLines(batext, con)
          
                ##Foreground rank representation
          
                figure = figure +1
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
                ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></A><center><b>Figure %s</b></center></tr></td></table>\n", "Foreground representation on the array", fpdf, fpng, figure)
                writeLines(ftext, con)
                legendlocal = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s:</b> False color representations of the spatial intensity (background and foreground) distributions of each arrays. The color scale is shown in the panel on the right. The color scale was chosen proportional to the ranks. These graphical representation permit to show problems during the experimentation such as fingerprints, artifactual gradient or dye specific failure for instance.</DIV>", figure)          
                writeLines(legendlocal, con)
              }

#################################
###Section 3 : Reproducibility###
#################################

            ##Scatter plots between replicates
            ##Have to find a better solution for the xlim and ylim.
            if("replicates" %in% names(phenoData(expressionset)))
              {
                section = section + 1
                sec3text = sprintf("<hr><h2>Section %s: Reproducibility</h2>", section)          
                writeLines(sec3text, con)

                figure = figure +1
                makePlot(con=con, name = paste(outfile, "_replicates", sep = ""),
                         w=10, h=10, fun = function() {
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
                         }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></td>\n", title="Scatter plots between replicates", fig = figure)
                legendrepli = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> represents the comparison of the log2 ratio of the replicate arrays. Each pairs of replicates are directly compared thanks to a scatter plot representation and a correlation coefficient is calculated. Scatterplot with a wide distribution of the dots and with a low correlation coefficient means a bad reproducibility of the experiment.</DIV>", figure)          
                writeLines(legendrepli, con)

              }
      
#################################################
###Section 4 : Homogeneity between experiments###
#################################################

            ##Boxplots
            section = section + 1
            sec4text = sprintf("<hr><h2>Section %s: Homogeneity between experiments</h2>", section)
            writeLines(sec4text, con)
            
            figure = figure + 1
            colours = brewer.pal(12, "Paired")
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

            makePlot(con=con, name = paste(outfile, "_boxplot", sep = ""),
                     w=8, h=8, fun = function() {
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
                     }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %sa</center></b></td>\n", title="Boxplots", fig = figure)
      
            ##Density if 1 group
            xlim = c(min(na.omit(dat)),max(na.omit(dat)))
            xlimr = c(min(na.omit(rc)),max(na.omit(rc)))
            xlimg = c(min(na.omit(gc)),max(na.omit(gc)))
      
            if(max(group) == 1)
              {
                makePlot(con=con, name = paste(outfile, "_density", sep = ""),
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
                           multi("ecdf",ldat,xlim,"","log(ratio)","")}, text="<tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b><Figure %sb</b></center></td></table>\n", title="Density plots", fig = figure)
              }
      
            ##Density if more than 1 group
            if(max(group) > 1)
              {
                dpng = paste(outfile,"_density.png", sep = "")
                dpdf = paste(outfile,"_density.pdf", sep = "")
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
                dtext = sprintf("<tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A></center></tr></td></table>\n", "Density plots", dpdf, dpng)
                writeLines(dtext, con)
              }
            legendhom = sprintf("<DIV STYLE=\"text-align:justify;\">The quality metrics in this section look at the distribution of the (raw, unnormalized) feature intensities for each array. <b>Figure %sa</b> shows density estimates (histograms), and <b>Figure %sb</b> presents boxplots of the same data. Arrays whose distributions are very different from the others should be considered for possible problems.</DIV>", figure, figure)          
            writeLines(legendhom, con)


########################################
###Section 5 : Array platform quality###
########################################
     
            ##GC
    
            if("GC" %in% rownames(featureData(expressionset)@varMetadata))
              {
                section = section + 1
                sec5text = sprintf("<hr><h2>Section %s: Array platform quality</h2>", section)    
                writeLines(sec5text, con)
            
                figure = figure + 1
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
          
                gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></tr></td></table>\n", "GC content effect ", gcpdf, "per array", " and ", gcopdf, "global", gcopdf, gcopng, figure)
                writeLines(gctext, con)
                legendgc = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> shows the distributions of the log2 intensities grouped by the percentage of cytosines (C) and guanines (G) among the nucleotides in each probe. Box plots (a), empirical cumulative distribution functions (ECDF, b) and kernel density estimates (c). Box and line colors in the three panels correspond to the same groups. Cytosine and guanine are able to form three hydrogen bonds, while adenine (A) and thymine (T) only form two, hence oligonucleotides with a higher proportion of C and G can form more stable hybridization bindings. This should result in higher intensities measured on the array, regardless of the abundance of target molecules.</DIV>",  figure)          
                writeLines(legendgc, con)
              }

            ##Mapping of probes
    
            if("gn" %in% rownames(featureData(expressionset)@varMetadata))
              {
                if(!"GC" %in% rownames(featureData(expressionset)@varMetadata))
                  {
                    section = section + 1
                    sec5text = sprintf("<hr><h2>Section %s: Array platform quality</h2>", section)          
                    writeLines(sec5text, con)
                  }
            
                figure = figure + 1

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

                gotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></tr></td></table>\n", "Gene mapping ", gpdf, "per array", " and ", gopdf, "global", gopdf, gopng, figure)
                writeLines(gotext, con)
                legendgo = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> shows the density distributions of the log2 ratios grouped by the mapping of the probes. Blue, density estimate of intensities of segments that overlap an annotated mRNA in RefSeq database. Gray, segments that do not overlap a Refseq mRNA.</DIV>",  figure)          
                writeLines(legendgo, con)


              }
            
##########################################
###Section 6 : Between array comparison###
##########################################
            ##Heatmap
            section = section + 1
            sec6text = sprintf("<hr><h2>Section %s: Between array comparison</h2>", section)
            writeLines(sec6text, con)
            figure = figure + 1
            colourRange = brewer.pal(9,"Greys")
            marg = ceiling(long/(sqrt(long)/1.5)+log10(numArrays))
      
            makePlot(con=con, name = paste(outfile, "_heatmap", sep = ""),
                     w=8, h=8, fun = function() {
                       heatmap(as.matrix(outM),
                               labRow = sN,
                               labCol = sN,
                               col = colourRange,
                               scale = "none",
                               main = "",
                               margins = c(marg,marg))
                     }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td>\n", title="Heatmap representation of the distance between experiments", fig = figure)
      
            m = matrix(pretty(outM,9),nrow=1,ncol=length(pretty(outM,9)))
            hlpng = paste(outfile, "_heatmaplegend.png", sep = "")  
            png(file= hlpng, width = 150, height = 450)
            image(m,xaxt="n",yaxt="n",ylab="Distance", col = colourRange, cex.lab = 0.8, mgp = c(1.5,1,0) )
            axis(2, label= as.list(pretty(outM,9)),at=seq(0,1,by=(1/(length(pretty(outM,9))-1))), cex.axis = 0.8, padj = 1)
            dev.off()
            hltext = sprintf("<td><IMG BORDER = \"0\" SRC=\"%s\"/></td></tr></table>\n",  hlpng)
            writeLines(hltext, con)
            legendheatmap = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> shows a false color heatmap of between arrays distances, computed as the MAD of the M-value for each pair of arrays. <br>dij = c.median|xmi-xmj|<br> Here, xmi is the normalized intensity value of the m-th probe on the i-th array, on the original data scale. c = 1.4826 is a constant factor that ensures consistency with the empirical variance for Normally distributed data (see manual page of the mad function in R). This plot can serve to detect outlier arrays. Consider the following decomposition of xmi: xmi = zm + Bmi + Emi, (1)where zm is the probe effect for probe m (the same across all arrays), Emi are i.i.d. random variables with mean zero and Bmi is such that for any array i, the majority of values Bmi are negligibly small (i. e. close to zero). Bmi represents differential expression effects. In this model, all values dij are (in expectation) the same, namely sqrt(2) times the standard deviation of Emi . Arrays whose distance matrix entries are way different give cause for suspicion. The dendrogram on this plot also can serve to check if, without any probe filtering, the experiments cluster accordingly to a biological meaning.</DIV>",  figure)          
            writeLines(legendheatmap, con)

      
            writeLines("</table>", con)
            closeHtmlPage(con)
    
          }##endofarrayQualityMetrics
          )

##plotting functions for ExpressionSet

aqm.expressionset = function(expressionset, outfile, log.transformed, numberofgraphs)
  {
   
    fn  = paste(outfile,"_QMreport",sep="")
    title = paste(outfile, " quality metrics report", sep="")
    titletext = sprintf("<hr><h1><center>%s</h1></center><table border = \"0\" cellspacing = 5 cellpadding = 2>", title)
    con = openHtmlPage(fn, title)
    writeLines(titletext, con)

    ##data preparation
    if(log.transformed == FALSE)
      dat = log2(exprs(expressionset))
    if(log.transformed == TRUE)
      dat = exprs(expressionset)
   
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
    
    outM = as.dist(dist2(na.omit(dat)))
    hc = hclust(outM)
    
    if(!missing(numberofgraphs))
      k = numberofgraphs
    if(missing(numberofgraphs) && numArrays > 50)
      k = ceiling(numArrays/50)
    if(missing(numberofgraphs) && numArrays <= 50)
      k = 1    
    group = cutree(hc, k = k)

##########################################
###Section 1 : Individual array quality###
##########################################
    
    ##MAplots
    ##function from affyQCReport
    section = 1
    figure = 1
    sec1text = sprintf("<hr><h2>Section %s: Individual array quality</h2>", section)
    writeLines(sec1text, con)
    
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
    
    matext1 = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></td><td><A HREF=\"%s\">%s</td></A></tr></td>\n", "MvA plots",  mapdf[1],  mapng[1], figure , mapdf[1], "MvA plot 1")
    writeLines(matext1, con)
    if(nfig >= 2)
      {
        for(i in 2:nfig)
          {
            matext2 = sprintf("<tr><td></td><td></td><td><A HREF=\"%s\">%s%s</A></CENTER><BR></tr></td>\n", mapdf[i], "MvA plot ", i)
            writeLines(matext2, con)
          }
      }
    writeLines("</table>", con)
    legendMA = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> represents MA plot for each array. <br> M A-plots are useful for pairwise comparisons between arrays. M and A are defined as :<br>
M = log2(I1) - log2(I2)<br>
A = 1/2 log2(I1)+log2(I2))<br>
where I1 and I2 are the vectors of normalized intensities of two channels, on the original data scale (i. e. not logarithm-transformed). Rather than comparing each array to every other array, here we compare each array to a single median  \"pseudo\"-array. Typically, we expect the mass of the distribution in an M A-plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A. Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).</DIV>", figure)          
    writeLines(legendMA, con)
                
#################################
###Section 2 : Spatial Effects###
#################################

#################################
###Section 3 : Reproducibility###
#################################

    ##Scatter plots between replicates
    ##Have to find a better solution for the xlim and ylim.
    if("replicates" %in% names(phenoData(expressionset)))
      {
        section = section + 1
        sec3text = sprintf("<hr><h2>Section %s: Reproducibility</h2>", section)          
        writeLines(sec3text, con)

        figure = figure +1
        makePlot(con=con, name = paste(outfile, "_replicates", sep = ""),
                 w=10, h=10, fun = function() {
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
                 }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></td>\n", title="Scatter plots between replicates", fig = figure)
                legendrepli = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> represents the comparison of the log2 intensities of the replicate arrays. Each pairs of replicates are directly compared thanks to a scatter plot representation and a correlation coefficient is calculated. Scatterplot with a wide distribution of the dots and with a low correlation coefficient means a bad reproducibility of the experiment.</DIV>", figure)          
                writeLines(legendrepli, con)
      }


#################################################
###Section 4 : Homogeneity between experiments###
#################################################

    ##Boxplots
    section = section + 1
    sec4text = sprintf("<hr><h2>Section %s: Homogeneity between experiments</h2>", section)
    writeLines(sec4text, con)
            
    figure = figure + 1
    makePlot(con=con, name = paste(outfile, "_boxplot", sep = ""),
             w=8, h=8, fun = function() {
               colours = brewer.pal(12, "Paired")
               par(cex.axis = 0.8, pty = "s", lheight =((1/log10(numArrays))*long), mai = c(((long/12)+0.2),0.4,0.2,0.2) , omi = c(0,0,0,0))
    boxplot(ldat, col = colours[2], las = 3, range = 0, names = sN)
             }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %sa</center></b></td>\n", title = "Boxplots", fig = figure)
    
    ##Density if 1 group
    xlim = c(min(na.omit(dat)),max(na.omit(dat)))
    if(max(group) == 1)
      {
        makePlot(con=con, name = paste(outfile, "_density", sep = ""),
                 w=10, h=10, fun = function() {
                   xlim = c(min(na.omit(dat)),max(na.omit(dat)))
                   nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(2.8,2.8),c(1.8,2), TRUE)
                   par(xaxt = "n", cex.axis = 0.8, mar = c(0,5,2,5))
                   multi("density",ldat,xlim,"Density","","")
                   par(xaxt = "s", cex.axis = 0.8, mar = c(4,5,0,5))
                   multi("ecdf",ldat,xlim,"ECDF","log(intensity)","")
                 }, text="<tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b><Figure %sb</b></center></td></table>\n", title="Density plots", fig = figure)
      }
                
    ##Density if more than 1 group
    if(max(group) > 1)
      {
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
        dtext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><b>Figure %sb</b></CENTER><BR></tr></td></table>\n", "Density plots", dpdf, dpng, figure)
        writeLines(dtext, con)
      }
    legendhom = sprintf("<DIV STYLE=\"text-align:justify;\">The quality metrics in this section look at the distribution of the (raw, unnormalized) feature intensities for each array. <b>Figure %sa</b> shows density estimates (histograms), and <b>Figure %sb</b> presents boxplots of the same data. Arrays whose distributions are very different from the others should be considered for possible problems.</DIV>", figure, figure)          
    writeLines(legendhom, con)



########################################
###Section 5 : Array platform quality###
########################################

        ##GC
        if("GC" %in% rownames(featureData(expressionset)@varMetadata))
          {
            section = section + 1
            sec5text = sprintf("<hr><h2>Section %s: Array platform quality</h2>", section)    
            writeLines(sec5text, con)
            figure = figure + 1
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
                
            gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></tr></td></table>\n", "GC content effect ", gcpdf, "per array", " and ", gcopdf, "global", gcopdf, gcopng, figure)
            writeLines(gctext, con)
            legendgc = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> shows the distributions of the log2 intensities grouped by the percentage of cytosines (C) and guanines (G) among the nucleotides in each probe. Box plots (a), empirical cumulative distribution functions (ECDF, b) and kernel density estimates (c). Box and line colors in the three panels correspond to the same groups. Cytosine and guanine are able to form three hydrogen bonds, while adenine (A) and thymine (T) only form two, hence oligonucleotides with a higher proportion of C and G can form more stable hybridization bindings. This should result in higher intensities measured on the array, regardless of the abundance of target molecules.</DIV>",  figure)          
            writeLines(legendgc, con)
          }

        ##Mapping of probes
        if("gn" %in% rownames(featureData(expressionset)@varMetadata))
          {
            if(!"GC" %in% rownames(featureData(expressionset)@varMetadata))
              {
                section = section + 1
                sec5text = sprintf("<hr><h2>Section %s: Array platform quality</h2>", section)          
                writeLines(sec5text, con)
              }
            
            figure = figure + 1
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
            
            gotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></tr></td></table>\n", "Gene mapping ", gpdf, "per array", " and ", gopdf, "global", gopdf, gopng, figure)
            writeLines(gotext, con)
            legendgo = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> shows the density distributions of the log2 ratios grouped by the mapping of the probes. Blue, density estimate of intensities of segments that overlap an annotated mRNA in RefSeq database. Gray, segments that do not overlap a Refseq mRNA.</DIV>",  figure)          
            writeLines(legendgo, con)
          }

    
##########################################
###Section 6 : Between array comparison###
##########################################
   
    ##Heatmap
    section = section + 1
    sec6text = sprintf("<hr><h2>Section %s: Between array comparison</h2>", section)
    writeLines(sec6text, con)
    figure = figure + 1
    colourRange = brewer.pal(9,"Greys")
    marg = ceiling(long/(sqrt(long)/1.5)+log10(numArrays))
    
    makePlot(con=con, name = paste(outfile, "_heatmap", sep = ""),
             w=8, h=8, fun = function() {
               heatmap(as.matrix(outM),
                       labRow = sN,
                       labCol = sN,
                       col = colourRange,
                       scale = "none",
                       main = "",
                       margins = c(marg,marg))
             }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td>\n", title="Heatmap representation of the distance between experiments", fig = figure)
      
    m = matrix(pretty(outM,9),nrow=1,ncol=length(pretty(outM,9)))
    hlpng = paste(outfile, "_heatmaplegend.png", sep = "")  
    png(file= hlpng, width = 150, height = 450)
    image(m,xaxt="n",yaxt="n",ylab="Distance", col = colourRange, cex.lab = 0.8, mgp = c(1.5,1,0) )
    axis(2, label= as.list(pretty(outM,9)),at=seq(0,1,by=(1/(length(pretty(outM,9))-1))), cex.axis = 0.8, padj = 1)
    dev.off()
    hltext = sprintf("<td><IMG BORDER = \"0\" SRC=\"%s\"/></CENTER><BR></td></tr></table>\n",  hlpng)
    writeLines(hltext, con)
    legendheatmap = sprintf("<DIV STYLE=\"text-align:justify;\"><b>Figure %s</b> shows a false color heatmap of between arrays distances, computed as the MAD of the M-value for each pair of arrays. <br>dij = c.median|xmi-xmj|<br> Here, xmi is the normalized intensity value of the m-th probe on the i-th array, on the original data scale. c = 1.4826 is a constant factor that ensures consistency with the empirical variance for Normally distributed data (see manual page of the mad function in R). This plot can serve to detect outlier arrays. Consider the following decomposition of xmi: xmi = zm + Bmi + Emi, (1)where zm is the probe effect for probe m (the same across all arrays), Emi are i.i.d. random variables with mean zero and Bmi is such that for any array i, the majority of values Bmi are negligibly small (i. e. close to zero). Bmi represents differential expression effects. In this model, all values dij are (in expectation) the same, namely sqrt(2) times the standard deviation of Emi . Arrays whose distance matrix entries are way different give cause for suspicion. The dendrogram on this plot also can serve to check if, without any probe filtering, the experiments cluster accordingly to a biological meaning.</DIV>",  figure)          
    writeLines(legendheatmap, con)

    l = list(numArrays,sN, section, figure, con)
    return(l) 
  }

##ExpressionSet
setMethod("arrayQualityMetrics",signature(expressionset="ExpressionSet"),
          function(expressionset, outfile, log.transformed, numberofgraphs)
          {
            l = aqm.expressionset(expressionset, outfile, numberofgraphs)
            con = l[[5]]
            writeLines("</table>", con)
            closeHtmlPage(con)

          })##end set method ExpressionSet



##AffyBatch

setMethod("arrayQualityMetrics",signature(expressionset="AffyBatch"),
          function(expressionset, outfile, log.transformed, numberofgraphs){
            
############################
###Section 7 : Affy plots###
############################
            
            l = aqm.expressionset(expressionset, outfile, numberofgraphs)
            numArrays = as.numeric(l[1])
            sN = l[[2]]
            section = as.numeric(l[[3]])
            figure = as.numeric(l[[4]])
            con = l[[5]]
            
            section = section + 1
            sec7text = sprintf("<hr><h2>Section %s: Affy plots</h2>", section)
            writeLines(sec7text, con)
            
            figure1 = figure + 1
            acol = sample(brewer.pal(8, "Dark2"), numArrays, replace = (8<numArrays))
            rnaDeg = AffyRNAdeg(expressionset)
            affypng1 = paste(outfile,"_RNAdeg.png", sep = "")
            affypdf1 = paste(outfile,"_RNAdeg.pdf", sep = "")
            png(file = affypng1)
            plotAffyRNAdeg(rnaDeg, cols = acol, lwd = 2)
            dev.copy(pdf, file = affypdf1)
            dev.off()
            dev.off()
                
            figure2 = figure1 + 1
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
            
            figure3 = figure2 + 1
            affypng3 = paste(outfile,"_NUSE.png", sep = "")
            affypdf3 = paste(outfile,"_NUSE.pdf", sep = "")
            png(file = affypng3)
            boxplot(dataPLM, ylim = c(0.95, 1.5), names = sN,
                    outline = FALSE, col = "lightblue", main = "NUSE", las = 2)
            dev.copy(pdf, file = affypdf3)
            dev.off()
            dev.off()
            affytext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</b></CENTER></tr></td></table>\n", "RNA degradation plot", affypdf1, affypng1, figure1, "RLE plot", affypdf2, affypng2, figure2, "NUSE plot", affypdf3, affypng3, figure3)
            writeLines(affytext, con)
            legendaffy = sprintf("<DIV STYLE=\"text-align:justify;\">In this section we present diagnostic plots based on tools provided in the affyPLM package.
In <b>Figure %s</b> a RNA digestion plot is computed. In this plot each array is represented by a single line. It is important to identify any array(s) that has a slope which is very different from the others. The indication is that the RNA used for that array has potentially been handled quite differently from the other arrays.
<b>Figure %s</b> is a Normalized Unscaled Standard Error (NUSE) plot. Low quality arrays are those that are significantly elevated or more spread out, relative to the other arrays. NUSE values are not comparable across data sets.
<b>Figure %s</b> is a Relative Log Expression (RLE) plot and an array that has problems will either have larger spread, or will not be centered at M = 0, or both.
</DIV>",  figure1, figure2, figure3)          
            writeLines(legendaffy, con)

            
            writeLines("</table>", con)
            closeHtmlPage(con)
          })##end set method AffyBatch

