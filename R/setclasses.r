setClass("aqmobj.rle", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.nuse", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.rnadeg", representation(plot="ANY", type="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.probesmap", representation(plot="ANY", type="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.qcs", representation(plot="ANY", type="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.pmmm", representation(plot="list", type="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.box", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.dens", representation(plot="ANY", type="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.heat", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.ma", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.spatial", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.spatialbg", representation(plot="ANY", type="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.msd", representation(plot="ANY", type="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.prepdata", representation(M = "matrix", A = "matrix", dat = "matrix", outM = "ANY", sN = "ANY", numArrays = "numeric", nchannels = "numeric", logtransformed = "logical", classori = "character"))

setClass("aqmobj.prepaffy", representation(dataPLM = "PLMset", sN = "character"))
