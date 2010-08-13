setClass("aqmobj.rle", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.nuse", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.rnadeg", representation(plot="ANY", section="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.probesmap", representation(plot="ANY", section="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.qcs", representation(plot="ANY", section="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.pmmm", representation(plot="list", section="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.box", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "list"))

setClass("aqmobj.dens", representation(plot="ANY", section="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.heat", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.pca", representation(plot="ANY", section="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.ma", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.spatial", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "character"))

setClass("aqmobj.spatialbg", representation(plot="ANY", section="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.msd", representation(plot="ANY", section="character", title="character", legend="character", shape = "character"))

setClass("aqmobj.prepdata", representation(rc = "matrix", gc = "matrix", rcb = "matrix", gcb = "matrix", M = "matrix", A = "matrix", dat = "matrix", outM = "ANY", sN = "ANY", numArrays = "numeric", nchannels = "numeric", logtransformed = "logical", classori = "character"))

setClass("aqmobj.prepaffy", representation(dataPLM = "PLMset", sN = "character"))
