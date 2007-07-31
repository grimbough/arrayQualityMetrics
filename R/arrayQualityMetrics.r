setGeneric("arrayQualityMetrics",
           function(expressionset,
                    outdir = getwd(),
                    force = FALSE,
                    do.logtransform = FALSE,
                    split.plots = FALSE)
           standardGeneric("arrayQualityMetrics")) 


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
  
  htmltext = sprintf(text, title, basename(outf[1]), basename(outf[2]), fig)
  
  resret = list(res,htmltext)
  
  return(resret)
}

##NChannelSet

setMethod("arrayQualityMetrics",signature(expressionset = "NChannelSet"),
          function(expressionset, outdir, force, do.logtransform, split.plots)
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

            olddir = getwd()
            setwd(outdir)

            arg = as.list(match.call(expand.dots = TRUE))
            title = paste(arg$expressionset, " quality metrics report", sep="")
            titletext = sprintf("<hr><h1><center>%s</h1></center><table border = \"0\" cellspacing = 5 cellpadding = 2 >", title)
            con = openHtmlPage("QMreport", title)
            writeLines(titletext, con)


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
                writeLines("<hr><h2><a name=\"EN\">Experiment Names</a></h2>", con)
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
      
            k = if(split.plots) split.plots else k = numArrays
            ##attribute randomly the experiments to different groups
            group = sample(rep((1:ceiling(numArrays/k)),k),numArrays)

##########################################
###Section 1 : Individual array quality###
##########################################
            
            section = 1
            figure = 1
            sec1text = sprintf("<hr><h2><a name=\"S1\">Section %s: Individual array quality</a></h2>", section)
            
           ##MAplots
            ##function from affyQCReport
            M = rc - gc
            A = 0.5*(rc + gc)
      
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
    
            matext1 = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></td><td><A HREF=\"%s\">%s</td></A></tr></td>\n", "MvA plots",  basename(mapdf[1]),  basename(mapng[1]), figure , basename(mapdf[1]), "MvA plot 1")
            legendMA = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> represents MA plot for each array. <br> M A-plots are useful for pairwise comparisons between arrays. M and A are defined as :<br>
M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>
A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>))<br>
where I<sub>2</sub> and I<sub>2</sub> are the vectors of normalized intensities of two channels, on the original data scale (i. e. not logarithm-transformed).Typically, we expect the mass of the distribution in an M A-plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A.
Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).</DIV>", figure)          
            
###############################
###Section 2 : Spatial plots###
###############################

            ##Background rank representation
      
            if("Row" %in% rownames(featureData(expressionset)@varMetadata) && "Column" %in% rownames(featureData(expressionset)@varMetadata))
              {
                section = section + 1
                sec2text = sprintf("<hr><h2><a name = \"S2\">Section %s: Spatial plots</a></h2>", section)

                r = featureData(expressionset)$Row
                c = featureData(expressionset)$Column

                maxc = max(as.numeric(c))
                maxr = max(as.numeric(r))
                colourRamp <- rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))

        
                if("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))
                  {
                    bapng = "background.png"
                    bapdf = "background.pdf"
                    pdf(file = bapdf)
                    nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = FALSE),
                                 c(1.2,1.2,1.2), c(2,2), TRUE)
                    for(a in 1:numArrays)
                      {
                        intr = cbind(r,c,rbg[,a])
                        intg = cbind(r,c,gbg[,a])
                  
                        re = matrix(NA,ncol=maxc,nrow=maxr)
                        g = matrix(NA,ncol=maxc,nrow=maxr)
                  
                        for(i in 1:nrow(intr))
                          {
                            re[intr[i,1],intr[i,2]] = intr[i,3]
                            g[intg[i,1],intg[i,2]] =  intg[i,3]
                          }
                  
                        mr = matrix(rank(re),ncol=maxc,nrow=maxr)
                        mg = matrix(rank(g),ncol=maxc,nrow=maxr)
                        if(maxr>maxc){
                          mr = t(mr)
                          mg = t(mg)
                        }
                        Imr = Image(mr, dim(mr))
                        Imrr = resize(normalize(Imr), w=nrow(mr)/3, h=ncol(mr)/3)
                        Img = Image(mg, dim(mg))
                        Imgr = resize(normalize(Img), w=nrow(mg)/3, h=ncol(mg)/3)
                        par(xaxt = "n", yaxt = "n",mar=c(1,1,2,1))
                        image(Imrr, col = colourRamp)
                        mtext(sN[a],side = 3,adj = 0.5, padj = -1 ,cex = 0.7)
                        mtext("Red Intensity",side = 2,adj = 0.5, padj = 4 ,cex = 0.7)
                        image(Imgr, col = colourRamp)
                        mtext("Green Intensity",side = 2,adj = 0.5, padj = 4 ,cex = 0.7)
                      }   
                    dev.off()

                    aR = ncol(Imrr)/nrow(Imrr)
                    w = 350
                    h = aR * w/3
                    png(file = bapng, width = w, height = h)
                    nf <- layout(matrix(c(1,2),1,2,byrow = FALSE),
                                 c(1.2,1.2), c(2,2), TRUE)
                    par(xaxt = "n", yaxt = "n",mar=c(1,1,1,1))
                    image(Imgr, col = colourRamp)
                    mtext("Red Intensity",side = 2,adj = 0.5, padj = 4 ,cex = 0.7)
                    image(Imgr, col = colourRamp)
                    mtext("Green Intensity",side = 2,adj = 0.5, padj = 4 ,cex = 0.7)
                    dev.off()
                    m = matrix(pretty(mr,9),nrow=1,ncol=length(pretty(mr,9)))
                    llbpng = "localisationlegendbackground.png"
                    png(file= llbpng, width = 150, height = 400)
                    image(m,xaxt="n",yaxt="n",ylab="Rank", col = colourRamp, cex.lab = 0.8, mgp = c(1.5,1,0) )
                    axis(2, label= as.list(pretty(mr,9)),at=seq(0,1,by=(1/(length(pretty(mr,9))-1))), cex.axis = 0.7, padj = 1)
                    dev.off()
                
                    batext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A></CENTER></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td></tr></table>\n", "Background representation on the array", basename(bapdf), basename(bapng), basename(llbpng))
                  }

                ##Foreground rank representation
          
                figure = figure +1
                fpng = "foreground.png"
                fpdf = "foreground.pdf"
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
                    
                    mrf = matrix(rank(rf),ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
                    mgf = matrix(rank(gf),ncol=max(as.numeric(c)),nrow=max(as.numeric(r)))
                    
                    if(maxr>maxc){
                      mrf = t(mrf)
                      mgf = t(mgf)
                    }
                    Imrf = Image(mrf, dim(mrf))
                    Imrfr = resize(normalize(Imrf), w=nrow(mrf)/3, h=ncol(mrf)/3)
                    Imgf = Image(mgf, dim(mgf))
                    Imgfr = resize(normalize(Imgf), w=nrow(mgf)/3, h=ncol(mgf)/3)
                    par(xaxt = "n", yaxt = "n",mar=c(1,1,2,1))
                    image(Imrfr, col = colourRamp)
                    mtext(sN[a],side = 3,adj = 0.5, padj = -1 ,cex = 0.7)
                    mtext("Red Intensity",side = 2,adj = 0.5, padj = 4 ,cex = 0.7)
                    image(Imgfr, col = colourRamp)
                    mtext("Green Intensity",side = 2,adj = 0.5, padj = 4 ,cex = 0.7)
                  }
                dev.off()
                
                aR = ncol(Imrfr)/nrow(Imrfr)
                w = 350
                h = aR * w/3
                png(file = fpng, width = w, height = h)
                nf <- layout(matrix(c(1,2),1,2,byrow = FALSE),
                             c(1.2,1.2), c(2,2), TRUE)
                par(xaxt = "n", yaxt = "n",mar=c(1,1,2,1))
                image(Imrfr, col = colourRamp)
                mtext("Red Intensity",side = 2,adj = 0.5, padj = 4 ,cex = 0.7)
                image(Imgfr, col = colourRamp)
                mtext("Green Intensity",side = 2,adj = 0.5, padj = 4 ,cex = 0.7)
                dev.off()
                
                m = matrix(pretty(mrf,9),nrow=1,ncol=length(pretty(mrf,9)))
                llfpng = "localisationlegendforeground.png"
                png(file= llfpng, width = 150, height = 400)
                image(m,xaxt="n",yaxt="n",ylab="Rank", col = colourRamp, cex.lab = 0.8, mgp = c(1.5,1,0) )
                axis(2, label= as.list(pretty(mrf,9)),at=seq(0,1,by=(1/(length(pretty(mrf,9))-1))), cex.axis = 0.7, padj = 1)
                dev.off()
                
                ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></A><center><b>Figure %s</b></center></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td></tr></table>\n", "Foreground representation on the array", basename(fpdf), basename(fpng), figure, basename(llfpng))

                legendlocal = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s:</b> False color representations of the spatial intensity (background and foreground) distributions of each arrays. The color scale is shown in the panel on the right. The color scale was chosen proportional to the ranks. These graphical representation permit to show problems during the experimentation such as fingerprints, artifactual gradient or dye specific failure for instance.</DIV>", figure)          
              }

#################################
###Section 3 : Reproducibility###
#################################
            
            ##Scatter plots between replicates
            ##To improve if we keep it.
            if(!TRUE)
              {
                if("Cy3" %in% rownames(featureData(expressionset)@varMetadata) && "Cy5" %in% rownames(featureData(expressionset)@varMetadata))
                  {
                    section = section + 1
                    sec3text = sprintf("<hr><h2><a name = \"S3\">Section %s: Reproducibility</a></h2>", section)          
                    figure = figure +1              

                    mplot1 = makePlot(con=con, name = "replicates",
                      w=10, h=10, fun = function() {
                        cy3cy5 = paste(expressionset@phenoData$Cy3,expressionset@phenoData$Cy5,sep="/")
                        lev = levels(as.factor(cy3cy5))
                        for(i in 1:length(lev))
                          {
                            subset = dat[, cy3cy5 == lev[i]]
                            comp = combn(ncol(subset),2)
                            layout(matrix(1:16,4,4,byrow = TRUE), rep(2,4), rep(2,4), TRUE)
                            par(mar = c(0,0,0,0))
                            for(j in 1:ncol(comp))
                              {
                                smoothScatter(subset[,comp[1,j]],subset[,comp[2,j]], xaxt = "n", xlim = xlim, yaxt = "n", ylim = xlim, xlab = "", ylab = "")
                                abline(a=0,b=1,col="red")
                              }
                          }
                      }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></td>\n", title="Scatter plots between replicates", fig = figure)
                    htmltext1 = mplot1[[2]]
                
                    legendrepli = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> represents the comparison of the log<sub>2</sub> ratio of the replicate arrays. Each pairs of replicates are directly compared thanks to a scatter plot representation and a correlation coefficient is calculated. Scatterplot with a wide distribution of the dots and with a low correlation coefficient means a bad reproducibility of the experiment.</DIV>", figure)          

                  }
              }
      
#################################################
###Section 4 : Homogeneity between experiments###
#################################################

            ##Boxplots
            section = section + 1
            sec4text = sprintf("<hr><h2><a name = \"S4\">Section %s: Homogeneity between experiments</a></h2>", section)
            
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

            mplot2 = makePlot(con=con, name = "boxplot",
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

            htmltext2 = mplot2[[2]]          
      
            ##Density if 1 group
            xlim = c(min(na.omit(dat)),max(na.omit(dat)))
            xlimr = c(min(na.omit(rc)),max(na.omit(rc)))
            xlimg = c(min(na.omit(gc)),max(na.omit(gc)))
      
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
                           multi("ecdf",ldat,xlim,"","log(ratio)","")}, text="<tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><center><br><b>Figure %sb</b></center></td></table>\n", title="Density plots", fig = figure)
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
                dtext = sprintf("<tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><center><br><b>Figure %sb</b></center></td></table>\n", "Density plots", basename(dpdf), basename(dpng), figure)
              }
            legendhom = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\">The quality metrics in this section look at the distribution of the feature intensities for each array. <b>Figure %sa</b> presents boxplots. On the left panel, the green boxes correspond to the log<sub>2</sub> intensities of the green channel and the red boxes to the log<sub>2</sub> intensities of the red channel. The right panel shows the boxplots of log<sub>2</sub>(ratio). <b>Figure %sb</b> shows density estimates (histograms) of the data. Arrays whose distributions are very different from the others should be considered for possible problems.</DIV>", figure, figure)          


########################################
###Section 5 : Array platform quality###
########################################
     
            ##GC
    
            if("GC" %in% rownames(featureData(expressionset)@varMetadata))
              {
                section = section + 1
                sec5text = sprintf("<hr><h2><a name = \"S5\">Section %s: Array platform quality</a></h2>", section)    
            
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
          
                gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></tr></td></table>\n", "GC content effect ", basename(gcpdf), "per array", " and ", basename(gcopdf), "global", basename(gcopdf), basename(gcopng), figure)
                legendgc = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the distributions of the log<sub>2</sub> intensities grouped by the percentage of cytosines (C) and guanines (G) among the nucleotides in each probe. Box plots (a), empirical cumulative distribution functions (ECDF, b) and kernel density estimates (c). Box and line colors in the three panels correspond to the same groups. Cytosine and guanine are able to form three hydrogen bonds, while adenine (A) and thymine (T) only form two, hence oligonucleotides with a higher proportion of C and G can form more stable hybridization bindings. This should result in higher intensities measured on the array, regardless of the abundance of target molecules.</DIV>",  figure)
              }

            ##Mapping of probes
    
            if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
              {
                if(!"GC" %in% rownames(featureData(expressionset)@varMetadata))
                  {
                    section = section + 1
                    sec5text = sprintf("<hr><h2><a name = \"S5\">Section %s: Array platform quality</a></h2>", section)          
                  }
            
                figure = figure + 1

                probemapping = expressionset@featureData$hasTarget
                facgene = as.vector(probemapping)          
                facgene[probemapping == "TRUE"] = 1
                facgene[probemapping == "FALSE"] = 0
                
                gpdf = "GenesMapping.pdf"
                gopng = "overall_GenesMapping.png"
                gopdf = "overall_GenesMapping.pdf"
          
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
                legend("topright", c("Mapped","Unmapped"),lty=1,lwd=2,col= c(cols[c(2,9)]), bty = "n")
                dev.copy(pdf, file = gopdf)
                dev.off()
                dev.off()

                gotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></tr></td></table>\n", "Gene mapping ", basename(gpdf), "per array", " and ", basename(gopdf), "global", basename(gopdf), basename(gopng), figure)
                legendgo = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the density distributions of the log<sub>2</sub> ratios grouped by the mapping of the probes. Blue, density estimate of intensities of probes annotated \"TRUE\" in the <b>\"hasTarget\"</b> slot. Gray, probes annotated \"FALSE\" in the <b>\"hasTarget\"</b> slot.</DIV>",  figure)          


              }
            
##########################################
###Section 6 : Between array comparison###
##########################################
            ##Heatmap
            section = section + 1
            sec6text = sprintf("<hr><h2><a name = \"S6\">Section %s: Between array comparison</a></h2>", section)
            figure = figure + 1
            colourRange = brewer.pal(9,"Greys")
            marg = ceiling(long/(sqrt(long)/1.5)+log10(numArrays))
      
            mplot4 = makePlot(con=con, name = "heatmap",
                     w=8, h=8, fun = function() {
                       heatmap(as.matrix(outM),
                               labRow = sN,
                               labCol = sN,
                               col = colourRange,
                               scale = "none",
                               main = "",
                               margins = c(marg,marg))
                     }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td>\n", title="Heatmap representation of the distance between experiments", fig = figure)

            htmltext4 = mplot4[[2]]
      
            m = matrix(pretty(outM,9),nrow=1,ncol=length(pretty(outM,9)))
            hlpng = "heatmaplegend.png"
            png(file= hlpng, width = 150, height = 450)
            image(m,xaxt="n",yaxt="n",ylab="Distance", col = colourRange, cex.lab = 0.8, mgp = c(1.5,1,0) )
            axis(2, label= as.list(pretty(outM,9)),at=seq(0,1,by=(1/(length(pretty(outM,9))-1))), cex.axis = 0.8, padj = 1)
            dev.off()
            hltext = sprintf("<td><IMG BORDER = \"0\" SRC=\"%s\"/></td></tr></table>\n",  basename(hlpng))
            legendheatmap = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows a false color heatmap of between arrays distances, computed as the MAD of the M-value for each pair of arrays. <br><center><i>d<sub>ij</sub> = c &bull median|x<sub>mi</sub>-x<sub>mj</sub>|</i></center><br><br> Here, <i>x<sub>mi</sub></i> is the normalized intensity value of the <i>m</i>-th probe on the <i>i</i>-th array, on the original data scale. <br><i>c = 1.4826</i> is a constant factor that ensures consistency with the empirical variance for Normally distributed data (see manual page of the mad function in R). This plot can serve to detect outlier arrays. <br>Consider the following decomposition of <i>x<sub>mi</sub>: x<sub>mi</sub> = z<sub>m</sub> + &beta<sub>mi</sub> + &epsilon<sub>mi</sub></i>, where <i>z<sub>m</sub></i> is the probe effect for probe <i>m</i> (the same across all arrays), <i>&epsilon<sub>mi</sub></i> are i.i.d. random variables with mean zero and <i>&beta<sub>mi</sub></i> is such that for any array <i>i</i>, the majority of values <i>&beta<sub>mi</sub></i> are negligibly small (i. e. close to zero). <i>&beta<sub>mi</sub></i> represents differential expression effects. In this model, all values <i>d<sub>ij</sub></i> are (in expectation) the same, namely <sqrt>2</sqrt> times the standard deviation of <i>&epsilon<sub>mi</i></sub> . Arrays whose distance matrix entries are way different give cause for suspicion. The dendrogram on this plot also can serve to check if, without any probe filtering, the experiments cluster accordingly to a biological meaning.</DIV></table>",  figure)
            
##########################################
###Section 7 : Variance Mean Dependency###
##########################################
            ##meanSdplot
            section = section + 1
            sec7text = sprintf("<hr><h2><a name = \"S7\">Section %s: Variance mean dependency</a></h2>", section)
            figure = figure + 1
      
            mplot5 = makePlot(con=con, name = "meanSd",
                     w=8, h=8, fun = function() {
                       meanSdPlot(cbind(rc, gc)) 
                     }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td></table>\n", title="Standard deviation versus rank of the mean", fig = figure)

            htmltext5 = mplot5[[2]]
      
            legendsdmean = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\">For each feature, this plot on <b>Figure %s</b> shows the empirical standard deviation of the normalized and glog-transformed data
on the <i>y</i>-axis versus the rank of the mean on the <i>x</i>-axis. The red dots, connected by lines, show the running median of the standard deviation. It should be approximately horizontal, that is, show no substantial trend.</DIV>",  figure)
            

#########################
### Writing the report###
#########################
            
            
            writeLines("<hr><h2>Index</h2><table border = \"0\" cellspacing = 5 cellpadding = 2>", con)
            
            writeLines("<tr><td><b><a href=\"#S1\">Individual array Quality</b></a></td></tr>", con)
            
            if("Row" %in% rownames(featureData(expressionset)@varMetadata) && "Column" %in% rownames(featureData(expressionset)@varMetadata))
              writeLines("<tr><td><b><a href=\"#S2\">Spatial plots</b></a></td></tr>", con)
            
            if(!TRUE)
              {
                if("replicates" %in% names(phenoData(expressionset)))
                  writeLines("<tr><td><b><a href=\"#S3\">Reproducibility</b></a></td></tr>", con) 
              }
            
            writeLines( "<tr><td><b><a href=\"#S4\">Homogeneity between experiments</b></a></td></tr>", con)
            
            if("GC" %in% rownames(featureData(expressionset)@varMetadata) || "hasTarget" %in% rownames(featureData(expressionset)@varMetadata))

              writeLines("<tr><td><b><a href=\"#S5\">Array platform quality</b></a></td></tr>", con)
            
            writeLines("<tr><td><b><a href=\"#S6\">Between array comparison</b></a></td></tr>" , con)
            writeLines("<tr><td><b><a href=\"#S7\">Variance mean dependency</b></a></td></tr>" , con)

            
            writeLines("</table><table border = \"0\" cellspacing = 5 cellpadding = 2>", con)

            writeLines(sec1text, con)
            writeLines(matext1, con)
            if(nfig >= 2)
              {
                for(i in 2:nfig)
                  {
                    matext2 = sprintf("<tr><td></td><td></td><td><A HREF=\"%s\">%s%s</A><BR></tr></td>\n", basename(mapdf[i]), "MvA plot ", i)
                    
                    writeLines(matext2, con)
                  }
              }
            writeLines("</table>", con)
            writeLines(legendMA, con)
            
            
            if("Row" %in% rownames(featureData(expressionset)@varMetadata) && "Column" %in% rownames(featureData(expressionset)@varMetadata))
              {
                
                writeLines(sec2text, con)
                  
                if("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))
                  writeLines(batext, con)
                
                writeLines(ftext, con)
                writeLines(legendlocal, con)
              }
            if(!TRUE)
              {
                if("replicates" %in% names(phenoData(expressionset)))
                  {
                    writeLines(sec3text, con)
                    writeLines(htmltext1, con)
                    writeLines(legendrepli, con)
                  }
              }

            writeLines(sec4text, con)
            writeLines(htmltext2, con)
            
            if(max(group) == 1)
              writeLines(htmltext3, con)
           
            if(max(group) > 1)
              writeLines(dtext, con)
              
            writeLines(legendhom, con)
            if("GC" %in% rownames(featureData(expressionset)@varMetadata))
              {

                writeLines(sec5text, con)
                writeLines(gctext, con)
                writeLines(legendgc, con)
              }

            if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
              {
                if(!"GC" %in% rownames(featureData(expressionset)@varMetadata))
                  {
                    writeLines(sec5text, con)
                  }
                writeLines(gotext, con)
                writeLines(legendgo, con)
              }
            
            writeLines(sec6text, con)
            writeLines(htmltext4, con)
            writeLines(hltext, con)            
            writeLines(legendheatmap, con)
            writeLines(sec7text, con)
            writeLines(htmltext5, con)
            writeLines(legendsdmean, con)            

            writeLines("</table>", con)
            closeHtmlPage(con)
            setwd(olddir)

    
          }##endofarrayQualityMetrics
          )

##plotting functions for ExpressionSet

aqm.expressionset = function(expressionset, outdir = getwd(), force = FALSE, do.logtransform = FALSE, split.plots = FALSE, arg)
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

    olddir = getwd()
    setwd(outdir)

    title = paste(arg$expressionset, " quality metrics report", sep="")
    titletext = sprintf("<hr><h1><center>%s</h1></center><table border = \"0\" cellspacing = 5 cellpadding = 2>", title)
    con = openHtmlPage("QMreport", title)
    writeLines(titletext, con)

   
    ##data preparation
    dat = if(do.logtransform) log2(exprs(expressionset)) else exprs(expressionset)
   
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
                col1 = "#e0e0f0"        #
                col2 = "#d0d0ff"        #
              }
            if(i %in% seq(2,length(sN),by = 2))
              {
                col1 = "#d0d0ff"        #
                col2 = "#e0e0f0"        #
              }
            stext = sprintf("<tr><td BGCOLOR=\"%s\"><b>%s</b></td><td BGCOLOR=\"%s\">%s</td></tr>", col1, sNt[i,2], col2, sNt[i,1])
            writeLines(stext, con)
          }
        writeLines("</table>", con)
      }
    
    outM = as.dist(dist2(na.omit(dat)))

    k = if(split.plots) split.plots else k = numArrays
    ##attribute randomly the experiments to different groups
    group = sample(rep((1:ceiling(numArrays/k)),k),numArrays)
   
##########################################
###Section 1 : Individual array quality###
##########################################
    
    ##MAplots
    ##function from affyQCReport
    section = 1
    figure = 1
    sec1text = sprintf("<hr><h2><a name = \"S1\">Section %s: Individual array quality</h2></a>", section)
    
    medArray = rowMedians(dat)
    M =  dat-medArray
    A = (dat+medArray)/2
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
    
    matext1 = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></td><td><A HREF=\"%s\">%s</td></A></tr></td>\n", "MvA plots",  basename(mapdf[1]),  basename(mapng[1]), figure , basename(mapdf[1]), "MvA plot 1")
    legendMA = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> represents MA plot for each array. <br> M A-plots are useful for pairwise comparisons between arrays. M and A are defined as :<br>
M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>
A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>))<br>
where I<sub>1</sub> and I<sub>2</sub> are the vectors of normalized intensities of two channels, on the original data scale (i. e. not logarithm-transformed). Rather than comparing each array to every other array, here we compare each array to a single median  \"pseudo\"-array. Typically, we expect the mass of the distribution in an M A-plot to be concentrated along the M = 0 axis, and there should be no trend in the mean of M as a function of A. Note that a bigger width of the plot of the M-distribution at the lower end of the A scale does not necessarily imply that the variance of the M-distribution is larger at the lower end of the A scale: the visual impression might simply be caused by the fact that there is more data at the lower end of the A scale. To visualize whether there is a trend in the variance of M as a function of A, consider plotting M versus rank(A).</DIV>", figure)          


#################################
###Section 3 : Reproducibility###
#################################

    ##Scatter plots between replicates
    ##To improve if we keep it.
    if(!TRUE)
      {
        if("replicates" %in% names(phenoData(expressionset)))
          {
            section = section + 1
            sec3text = sprintf("<hr><h2><a name = \"S3\">Section %s: Reproducibility</h2></a>", section)          
        
            figure = figure +1
            mplot1 = makePlot(con=con, name = "replicates",
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
            htmltext1 = mplot1[[2]]
            legendrepli = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> represents the comparison of the log<sub>2</sub> intensities of the replicate arrays. Each pairs of replicates are directly compared thanks to a scatter plot representation and a correlation coefficient is calculated. Scatterplot with a wide distribution of the dots and with a low correlation coefficient means a bad reproducibility of the experiment.</DIV>", figure)          
          }
      }

#################################################
###Section 4 : Homogeneity between experiments###
#################################################

    ##Boxplots
    section = section + 1
    sec4text = sprintf("<hr><h2><a name = \"S4\">Section %s: Homogeneity between experiments</h2></a>", section)
            
    figure = figure + 1
    mplot2 = makePlot(con=con, name = "boxplot",
      w=8, h=8, fun = function() {
        colours = brewer.pal(12, "Paired")
        par(cex.axis = 0.8, pty = "s", lheight =((1/log10(numArrays))*long), mai = c(((long/12)+0.2),0.4,0.2,0.2) , omi = c(0,0,0,0))
        boxplot(ldat, col = colours[2], las = 3, range = 0, names = sN)
      }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %sa</center></b></td>\n", title = "Boxplots", fig = figure)
    htmltext2 = mplot2[[2]]
    
    ##Density if 1 group
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
          }, text="<tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><center><br><b>Figure %sb</b></center></td></table>\n", title="Density plots", fig = figure)
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
        dtext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><center><br><b>Figure %sb</b></CENTER><BR></tr></td></table>\n", "Density plots", basename(dpdf), basename(dpng), figure)
      }
    legendhom = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\">The quality metrics in this section look at the distribution of the feature intensities for each array. <b>Figure %sa</b> presents boxplots and <b>Figure %sb</b> shows density estimates (histograms) of the data. Arrays whose distributions are very different from the others should be considered for possible problems.</DIV>", figure, figure)          




########################################
###Section 5 : Array platform quality###
########################################

    ##GC
    if("GC" %in% rownames(featureData(expressionset)@varMetadata))
      {
        section = section + 1
        sec5text = sprintf("<hr><h2><a name = \"S5\">Section %s: Array platform quality</h2></a>", section)    
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
                
        gctext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></tr></td></table>\n", "GC content effect ", basename(gcpdf), "per array", " and ", basename(gcopdf), "global", basename(gcopdf), basename(gcopng), figure)
        legendgc = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the distributions of the log<sub>2</sub> intensities grouped by the percentage of cytosines (C) and guanines (G) among the nucleotides in each probe. Box plots (a), empirical cumulative distribution functions (ECDF, b) and kernel density estimates (c). Box and line colors in the three panels correspond to the same groups. Cytosine and guanine are able to form three hydrogen bonds, while adenine (A) and thymine (T) only form two, hence oligonucleotides with a higher proportion of C and G can form more stable hybridization bindings. This should result in higher intensities measured on the array, regardless of the abundance of target molecules.</DIV>",  figure)          
      }

    ##Mapping of probes
    if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
      {
        if(!"GC" %in% rownames(featureData(expressionset)@varMetadata))
          {
            section = section + 1
            sec5text = sprintf("<hr><h2><a name = \"S5\">Section %s: Array platform quality</h2></a>", section)          
          }
            
        figure = figure + 1

        probemapping = expressionset@featureData$hasTarget
        facgene = as.vector(probemapping)          
        facgene[probemapping == "TRUE"] = 1
        facgene[probemapping == "FALSE"] = 0

                    
        gpdf =  "GenesMapping.pdf"
        gopng = "overall_GenesMapping.png"
        gopdf = "overall_GenesMapping.pdf"
        
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
        legend("topright", c("Mapped","Unmapped"),lty=1,lwd=2,col= c(cols[c(2,9)]), bty = "n")
        dev.copy(pdf, file = gopdf)
        dev.off()
        dev.off()
            
        gotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></tr></td></table>\n", "Gene mapping ", basename(gpdf), "per array", " and ", basename(gopdf), "global", basename(gopdf), basename(gopng), figure)
        legendgo = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the density distributions of the log<sub>2</sub> ratios grouped by the mapping of the probes. Blue, density estimate of intensities of probes annotated \"TRUE\" in the <b>\"hasTarget\"</b> slot. Gray, probes annotated \"FALSE\" in the <b>\"hasTarget\"</b> slot.</DIV>", figure)          
      }

    
##########################################
###Section 6 : Between array comparison###
##########################################
   
    ##Heatmap
    section = section + 1
    sec6text = sprintf("<hr><h2><a name = \"S6\">Section %s: Between array comparison</h2></a>", section)
    figure = figure + 1
    colourRange = brewer.pal(9,"Greys")
    marg = ceiling(long/(sqrt(long)/1.5)+log10(numArrays))
    
    mplot4 = makePlot(con=con, name = "heatmap",
      w=8, h=8, fun = function() {
        heatmap(as.matrix(outM),
                labRow = sN,
                labCol = sN,
                col = colourRange,
                scale = "none",
                main = "",
                margins = c(marg,marg))
      }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td>\n", title="Heatmap representation of the distance between experiments", fig = figure)
    htmltext4 = mplot4[[2]]
      
    m = matrix(pretty(outM,9),nrow=1,ncol=length(pretty(outM,9)))
    hlpng = "heatmaplegend.png"
    png(file= hlpng, width = 150, height = 450)
    image(m,xaxt="n",yaxt="n",ylab="Distance", col = colourRange, cex.lab = 0.8, mgp = c(1.5,1,0) )
    axis(2, label= as.list(pretty(outM,9)),at=seq(0,1,by=(1/(length(pretty(outM,9))-1))), cex.axis = 0.8, padj = 1)
    dev.off()
    hltext = sprintf("<td><IMG BORDER = \"0\" SRC=\"%s\"/></CENTER><BR></td></tr></table>\n",  basename(hlpng))
                legendheatmap = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows a false color heatmap of between arrays distances, computed as the MAD of the M-value for each pair of arrays. <br><center><i>d<sub>ij</sub> = c &bull median|x<sub>mi</sub>-x<sub>mj</sub>|</i></center><br><br> Here, <i>x<sub>mi</sub></i> is the normalized intensity value of the <i>m</i>-th probe on the <i>i</i>-th array, on the original data scale. <br><i>c = 1.4826</i> is a constant factor that ensures consistency with the empirical variance for Normally distributed data (see manual page of the mad function in R). This plot can serve to detect outlier arrays. <br>Consider the following decomposition of <i>x<sub>mi</sub>: x<sub>mi</sub> = z<sub>m</sub> + &beta<sub>mi</sub> + &epsilon<sub>mi</sub></i>, where <i>z<sub>m</sub></i> is the probe effect for probe <i>m</i> (the same across all arrays), <i>&epsilon<sub>mi</sub></i> are i.i.d. random variables with mean zero and <i>&beta<sub>mi</sub></i> is such that for any array <i>i</i>, the majority of values <i>&beta<sub>mi</sub></i> are negligibly small (i. e. close to zero). <i>&beta<sub>mi</sub></i> represents differential expression effects. In this model, all values <i>d<sub>ij</sub></i> are (in expectation) the same, namely <sqrt>2</sqrt> times the standard deviation of <i>&epsilon<sub>mi</i></sub> . Arrays whose distance matrix entries are way different give cause for suspicion. The dendrogram on this plot also can serve to check if, without any probe filtering, the experiments cluster accordingly to a biological meaning.</DIV></table>",  figure)

##########################################
###Section 7 : Variance Mean Dependency###
##########################################
    ##meanSdplot
    section = section + 1
    sec7text = sprintf("<hr><h2><a name = \"S7\">Section %s: Variance mean dependency</a></h2>", section)
    figure = figure + 1
      
    mplot5 = makePlot(con=con, name = "meanSd",
      w=8, h=8, fun = function() {
        meanSdPlot(dat) 
      }, text="<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><BR><b>Figure %s</center></b></td></table>\n", title="Standard deviation versus rank of the mean", fig = figure)

    htmltext5 = mplot5[[2]]
    
    legendsdmean = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\">For each feature, this plot on <b>Figure %s</b> shows the empirical standard deviation of the normalized and glog-transformed data
on the <i>y</i>-axis versus the rank of the mean on the <i>x</i>-axis. The red dots, connected by lines, show the running median of the standard deviation. It should be approximately horizontal, that is, show no substantial trend.</DIV>",  figure)
    

#########################
### Writing the report###
#########################
            
            
    writeLines("<hr><h2>Index</h2><table border = \"0\" cellspacing = 5 cellpadding = 2>", con)
            
    writeLines("<tr><td><b><a href=\"#S1\">Individual array Quality</b></a></td></tr>", con)
               
    if(!TRUE)
      {
        if("replicates" %in% names(phenoData(expressionset)))
          writeLines("<tr><td><b><a href=\"#S3\">Reproducibility</b></a></td></tr>", con)
      }
    
    writeLines( "<tr><td><b><a href=\"#S4\">Homogeneity between experiments</b></a></td></tr>", con)
    
    if("GC" %in% rownames(featureData(expressionset)@varMetadata) || "hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
      
      writeLines("<tr><td><b><a href=\"#S5\">Array platform quality</b></a></td></tr>", con)
    
    writeLines("<tr><td><b><a href=\"#S6\">Between array comparison</b></a></td></tr>" , con)
    writeLines("<tr><td><b><a href=\"#S7\">Variance mean dependency</b></a></td></tr>" , con)

    if(is(expressionset, "AffyBatch"))
      {
        writeLines("<tr><td><b><a href=\"#S2\">Spatial plots</b></a></td></tr><tr><td><b><a href=\"#S8\">Affymetrix specific plots</b></a></td></tr>" , con)
      }
    
    writeLines("</table><table border = \"0\" cellspacing = 5 cellpadding = 2>", con)
    
    writeLines(sec1text, con)
    writeLines(matext1, con)
    if(nfig >= 2)
      {
        for(i in 2:nfig)
          {
            matext2 = sprintf("<tr><td></td><td></td><td><A HREF=\"%s\">%s%s</A><BR></tr></td>\n", basename(mapdf[i]), "MvA plot ", i)
            
            writeLines(matext2, con)
          }
      }
    writeLines("</table>", con)
    writeLines(legendMA, con)    
    
    if(!TRUE)
      {
          if("replicates" %in% names(phenoData(expressionset)))
          {
            writeLines(sec3text, con)
            writeLines(htmltext1, con)
            writeLines(legendrepli, con)
          }
      }
    
    writeLines(sec4text, con)
    writeLines(htmltext2, con)
    if(max(group) == 1)
      writeLines(htmltext3, con)
    
    if(max(group) > 1)
      writeLines(dtext, con)
    
    writeLines(legendhom, con)
    if("GC" %in% rownames(featureData(expressionset)@varMetadata))
      {
        
        writeLines(sec5text, con)
        writeLines(gctext, con)
        writeLines(legendgc, con)
      }

    if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
      {
        if(!"GC" %in% rownames(featureData(expressionset)@varMetadata))
          {
            writeLines(sec5text, con)
          }
        writeLines(gotext, con)
        writeLines(legendgo, con)
      }
    writeLines(sec6text, con)
    writeLines(htmltext4, con)
    writeLines(hltext, con)            
    writeLines(legendheatmap, con)
    writeLines(sec7text, con)
    writeLines(htmltext5, con)
    writeLines(legendsdmean, con)            

    l = list(numArrays,sN, section, figure, con, dat, olddir)

    return(l) 
  }

##ExpressionSet
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



##AffyBatch

setMethod("arrayQualityMetrics",signature(expressionset="AffyBatch"),
          function(expressionset, outdir, force, do.logtransform,
          split.plots){
            
            arg = as.list(match.call(expand.dots = TRUE))

            l = aqm.expressionset(expressionset, outdir, force, do.logtransform, split.plots, arg)
            numArrays = as.numeric(l[1])
            sN = l[[2]]
            section = as.numeric(l[[3]])
            figure = as.numeric(l[[4]])
            con = l[[5]]
            dat = l[[6]]
            olddir = l[[7]]

                
###############################
###Section 2 : Spatial Plots###
###############################
  
         
            figure = figure + 1
            section = section + 1
            sec2text = sprintf("<hr><h2><a name = \"S2\">Section %s: Spatial plots</h2></a>", section)
            writeLines(sec2text, con)

            fpng = "foreground.png"
            fpdf = "foreground.pdf"
            colourRamp <- rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))
            pdf(file = fpdf)
            nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow = TRUE),
                         c(1.2,1.2,1.2), c(2,2), TRUE)
            for(a in 1:numArrays)
              {
                par(xaxt = "n", yaxt = "n",mar=c(1,1,1,1))
                rfi = rank(dat[,a])
                mrfi = matrix(rfi,ncol=ncol(expressionset),nrow=nrow(expressionset),byrow=T)
                Imrfi = Image(mrfi, dim(mrfi))
                mrir = resize(normalize(Imrfi), w=nrow(Imrfi)/3, h=ncol(Imrfi)/3)
                image(mrir, col = colourRamp)
                mtext(sN[a], side = 3, adj = 0.5, padj = 4 , cex = 0.8)
                
              }
            dev.off()

            aR = ncol(Imrfi)/nrow(Imrfi)
            w = 350
            h = aR * w
            png(file = fpng, width = w, height = h)
            par(xaxt = "n", yaxt = "n",mar=c(1,1,1,1))
            image(mrir, col = colourRamp)
            dev.off()
            
            m = matrix(pretty(mrir,9),nrow=1,ncol=length(pretty(mrir,9)))
            llpng = "localisationlegend.png" 
            png(file= llpng, width = 150, height = 400)
            image(m,xaxt="n",yaxt="n",ylab="Rank", col = colourRamp, cex.lab = 0.8, mgp = c(1.5,1,0) )
            axis(2, label= as.list(pretty(mrir,9)),at=seq(0,1,by=(1/(length(pretty(mrir,9))-1))), cex.axis = 0.7, padj = 1)
            dev.off()
            
            ftext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG border = \"0\" SRC=\"%s\"/></A><center><b>Figure %s</b></center></td><td><IMG BORDER = \"0\" SRC=\"%s\"/></td></tr></table>\n", "Intensity representation on the array", basename(fpdf), basename(fpng), figure, basename(llpng))
            
            writeLines(ftext, con)
            legendlocal = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s:</b> False color representations of the spatial intensity distributions of each arrays. The color scale is shown in the panel on the right. The color scale was chosen proportional to the ranks. These graphical representation permit to show problems during the experimentation such as fingerprints, artifactual gradient or dye specific failure for instance.</DIV>", figure)          
            writeLines(legendlocal, con)



############################
###Section 8 : Affy plots###
############################

            
            section = section + 1
            sec8text = sprintf("<hr><h2><a name = \"S8\">Section %s: Affymetrix specific plots</h2></a>", section)
            writeLines(sec8text, con)
            
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
            Mbox(dataPLM, ylim = c(-1, 1), names = sN, col = "lightblue",
                 whisklty = 0, staplelty = 0, main = "RLE", las = 3)
            dev.copy(pdf, file = affypdf2)
            dev.off()
            dev.off()
            
            figure3 = figure2 + 1
            affypng3 = "NUSE.png"
            affypdf3 = "NUSE.pdf"
            png(file = affypng3)
            boxplot(dataPLM, ylim = c(0.95, 1.5), names = sN,
                    outline = FALSE, col = "lightblue", main = "NUSE", las = 2)
            dev.copy(pdf, file = affypdf3)
            dev.off()
            dev.off()
            
            figure4 = figure3 + 1
            affypng4 = "qc.png"
            affypdf4 = "qc.pdf"
            png(file = affypng4)
            plot(qc(expressionset))
            dev.copy(pdf, file = affypdf4)
            dev.off()
            dev.off()

            affytext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><center><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><center><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><center><BR><b>Figure %s</b></CENTER></tr></td><tr><td><b>%s</b></td><td><A HREF=\"%s\"><IMG BORDER = \"0\" SRC=\"%s\"/></A><center><BR><b>Figure %s</b></CENTER></tr></td></table>\n", "RNA degradation plot", basename(affypdf1), basename(affypng1), figure1, "RLE plot", basename(affypdf2), basename(affypng2), figure2, "NUSE plot", basename(affypdf3), basename(affypng3), figure3, "Diagnostic plot recommended by Affy", basename(affypdf4), basename(affypng4), figure4)
            writeLines(affytext, con)
            legendaffy = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\">In this section we present diagnostic plots based on tools provided in the affyPLM package. In <b>Figure %s</b> a RNA digestion plot is computed. In this plot each array is represented by a single line. It is important to identify any array(s) that has a slope which is very different from the others. The indication is that the RNA used for that array has potentially been handled quite differently from the other arrays. <b>Figure %s</b> is a Normalized Unscaled Standard Error (NUSE) plot. Low quality arrays are those that are significantly elevated or more spread out, relative to the other arrays. NUSE values are not comparable across data sets. <b>Figure %s</b> is a Relative Log Expression (RLE) plot and an array that has problems will either have larger spread, or will not be centered at M = 0, or both. <b>Figure %s</b> represents the diagnostic plot recommended by Affymetrix. It is fully describe in the QCandsimpleaffy.pdf vignette of the package simpleaffy. Any metrics (circles and triangles) that is shown in red is out of the manufacturer's specific boundaries and suggests a potential problem.</DIV>", figure1, figure2, figure3, figure4)            
            writeLines(legendaffy, con)
            
            ##PM.MM
            figure5 = figure4 + 1
            
            cols = brewer.pal(9, "Set1")
            xlim = c(min(na.omit(dat)),max(na.omit(dat)))

            
            pmpdf =  "PM.MM.pdf"         
            pmopng = "overall_PM.MM.png"
            pmopdf = "overall_PM.MM.pdf"
        
            pdf(pmpdf)
            nf <- layout(matrix(1:16,4,4,byrow=TRUE), c(2,1.8,1.8,1.8), c(1.8,1.8,1.8,2), FALSE)
            for(a in 1:numArrays)
              {
                if(a %in% c(seq(14,numArrays,by=16),seq(15,numArrays,by=16),seq(16,numArrays,by=16)))
                  {
                    par(mar = c(2,0,0,0))
                    pmmm(expressionset[,a],xlim,"","","", yaxt = "n", ylim = c(0,1))
                    legend("topright",legend=sN[a], bty = "n", cex = 0.6)
                  }
                if(a %in%  c(seq(1,numArrays,by=16),seq(5,numArrays,by=16),seq(9,numArrays,by=16)))
                  {
                    par(mar = c(0,2,0,0))
                    pmmm(expressionset[,a],xlim,"","","", xaxt = "n", ylim = c(0,1))
                    legend("topright",legend=sN[a], bty = "n", cex = 0.6)
                  }
                if(a %in% seq(13,numArrays,by=16))
                  {
                    par(mar = c(2,2,0,0))
                    pmmm(expressionset[,a],xlim,"","","", ylim = c(0,1))
                    legend("topright",legend=sN[a], bty = "n", cex = 0.6)
                  }
                if(!(a %in%  seq(13,numArrays,by=16)) && !(a %in%  c(seq(1,numArrays,by=16),seq(5,numArrays,by=16),seq(9,numArrays,by=16))) && !(a %in% c(seq(14,numArrays,by=16),seq(15,numArrays,by=16),seq(16,numArrays,by=16))))
                  {
                    par(mar = c(0,0,0,0))
                    pmmm(expressionset[,a],xlim,"","","", xaxt = "n", yaxt = "n",ylim = c(0,1))
                    legend("topright",legend=sN[a], bty = "n", cex = 0.6)
                  }
              }
            dev.off()
        
            png(file = pmopng)
            pmmm(expressionset,xlim,"","","")
            legend("topright", c("PM","MM"),lty=1,lwd=2,col= c(cols[c(2,9)]), bty = "n")
            dev.copy(pdf, file = pmopdf)
            dev.off()
            dev.off()
            
            pmotext = sprintf("<table cellspacing = 5 cellpadding = 2><tr><td><b>%s<A HREF=\"%s\">%s</A>%s<A HREF=\"%s\">%s</A></b><td><A HREF=\"%s\"><center><IMG BORDER = \"0\" SRC=\"%s\"/></A><br><b>Figure %s</b></center></tr></td></table>\n", "Perfect matchs and mismatchs ", basename(pmpdf), "per array", " and ", basename(pmopdf), "global", basename(pmopdf), basename(pmopng), figure5)
            legendpmo = sprintf("<DIV style=\"font-size: 13; font-family: Lucida Grande; text-align:justify\"><b>Figure %s</b> shows the density distributions of the log<sub>2</sub> intensities grouped by the matching of the probes. Blue, density estimate of intensities of perfect match probes and gray the mismatch probes.</DIV>",  figure5)
            writeLines(pmotext, con)
            writeLines(legendpmo, con)

            
            writeLines("</table>", con)
            closeHtmlPage(con)
            setwd(olddir)

          })##end set method AffyBatch

