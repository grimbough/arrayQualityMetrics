## x : featureData(NChannelSet)
## gal.file : name of the file .gal to use
## nBlocks : number of blocks on the array
## skip : number of header lines to skip in the gal.file
## require(Biobase)

addXYfromGAL = function(x, gal.file, nBlocks, skip, ...)
  {
    galtype = readLines(gal.file,n=(skip+1))
    sep = if(length(grep("\t",galtype[skip+1])) == 1) "\t" else " "
    
    gal = read.table(gal.file, skip = skip, sep = sep, nrows = nBlocks, quote="", as.is = T, ...)
    
    if(length(grep(",",galtype[skip+1])) == 1)
      gal = sapply(1:ncol(gal),function(i) gsub(",","",gal[,i]))

    if(length(grep("=$",gal[,1])) == 0)
      gal = cbind(matrix(unlist(strsplit(gal[,1],"=")),ncol=2,byrow=T),gal[,2:ncol(gal)])
    
    block = seq_len(nBlocks)

    Xfac = seq_len(length(levels(as.factor(gal[,2]))))
    Yfac = seq_len(length(levels(as.factor(gal[,3]))))
    bcoord = cbind(block,Xfac,Yfac)

    r = x$Row
    c = x$Column
    b = x$Block

    bc = bcoord[b[],]
    absr = as.numeric(bc[,2])*r
    absc = as.numeric(bc[,3])*c
       
    x$X = absr
    x$Y = absc
    return(x)
  }

