setClass("aqmobj.rle", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY"))

setClass("aqmobj.nuse", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY"))

setClass("aqmobj.rnadeg", representation(plot="ANY", type="character", title="character", legend="character"))

setClass("aqmobj.probesmap", representation(plot="ANY", type="character", title="character", legend="character"))

setClass("aqmobj.qcs", representation(plot="ANY", type="character", title="character", legend="character"))

setClass("aqmobj.pmmm", representation(plot="list", type="character", title="character", legend="character"))

setClass("aqmobj.box", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY"))

setClass("aqmobj.dens", representation(plot="ANY", type="character", title="character", legend="character"))

setClass("aqmobj.heat", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY"))

setClass("aqmobj.ma", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY"))

setClass("aqmobj.spatial", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY"))

setClass("aqmobj.spatialbg", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY"))

setClass("aqmobj.msd", representation(plot="ANY", type="character", title="character", legend="character"))

setClass("aqmobj.prepdata", representation(M = "matrix", A = "matrix", dat = "matrix", outM = "ANY", sN = "ANY", numArrays = "numeric", nchannels = "numeric", logtransformed = "logical", classori = "character"))

setClass("aqmobj.prepaffy", representation(dataPLM = "PLMset", sN = "character"))
