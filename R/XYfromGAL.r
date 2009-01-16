## x : featureData(NChannelSet)
## gal.file : name of the file .gal to use
## nBlocks : number of blocks on the array
## skip : number of header lines to skip in the gal.file
## require(Biobase)
makeabscoord = function(a) {
  sapply(a, function(x) order(unique(a))[unique(a)==x])
}

makecoordblock = function(a, coo, bcoord) {
  sapply(seq_len(nrow(coo)), function(i) coo[i,a]*bcoord[coo[i,1],a])
}

addXYfromGAL = function(x, gal.file, nBlocks, skip, ...)
  {
    galtype = readLines(gal.file,n=(skip+1))
    sep = if(length(grep("\t",galtype[skip+1])) == 1) "\t" else " "
    
    gal = read.table(gal.file, skip = skip, sep = sep, nrows = nBlocks, quote="", as.is = T, ...)
    
    if(length(grep(",",galtype[skip+1])) == 1)
      gal = sapply(seq_len(ncol(gal)),function(i) gsub(",","",gal[,i]))

    if(length(grep("=$",gal[,1])) == 0)
      gal = cbind(matrix(unlist(strsplit(gal[,1],"=")),ncol=2,byrow=T),gal[,2:ncol(gal)])
    
    block = seq_len(nBlocks)

    Xfac = makeabscoord(as.numeric(gal[,2]))
    Yfac = makeabscoord(as.numeric(gal[,3]))
    bcoord = cbind(block,Xfac,Yfac)

    coo = cbind(x$Block, x$Column, x$Row)

    x$X = makecoordblock(2, coo, bcoord)
    x$Y = makecoordblock(3, coo, bcoord)
    return(x)
  }

