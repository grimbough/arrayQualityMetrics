setClass("aqmobj.rle", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "list"))

setClass("aqmobj.nuse", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "list"))

setClass("aqmobj.rnadeg", representation(plot="ANY", section="character", title="character", legend="character", shape = "list"))

setClass("aqmobj.probesmap", representation(plot="ANY", section="character", title="character", legend="character", shape = "list"))

setClass("aqmobj.qcs", representation(plot="ANY", section="character", title="character", legend="character", shape = "list"))

setClass("aqmobj.pmmm", representation(plot="list", section="character", title="character", legend="character", shape = "list"))

setClass("aqmobj.box", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "list"))

setClass("aqmobj.dens", representation(plot="ANY", section="character", title="character", legend="character", shape = "list"))

setClass("aqmobj.heat", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "list"))

setClass("aqmobj.pca", representation(plot="ANY", section="character", title="character", legend="character", shape = "list"))

setClass("aqmobj.ma", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "list"))

setClass("aqmobj.spatial", representation(plot="ANY", section="character", title="character", legend="character", scores="numeric", outliers="ANY", shape = "list"))

setClass("aqmobj.spatialbg", representation(plot="ANY", section="character", title="character", legend="character", shape = "list"))

setClass("aqmobj.msd", representation(plot="ANY", section="character", title="character", legend="character", shape = "list"))

setClass("aqmobj.prepdata", representation(rc = "matrix", gc = "matrix", rcb = "matrix", gcb = "matrix", M = "matrix", A = "matrix", dat = "matrix", outM = "ANY", sN = "ANY", numArrays = "numeric", nchannels = "numeric", logtransformed = "logical", classori = "character"))

setClass("aqmobj.prepaffy", representation(dataPLM = "PLMset"))
