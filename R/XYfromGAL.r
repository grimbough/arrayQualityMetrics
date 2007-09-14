## x : featureData(NChannelSet)
## gal.file : name of the file .gal to use
## nBlocks : number of blocks on the array
## skip : number of header lines to skip in the gal.file
## require(Biobase)

addXYfromGAL = function(x, gal.file, nBlocks, skip, ...)
  {
    gal = read.table(gal.file, skip = skip, sep = "\t", nrows = nBlocks, quote="", as.is = T, ...)
    block = gsub("\"Block", "", gal[,1])
    block = gsub("=", "", block)

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

