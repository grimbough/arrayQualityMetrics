## Suppose the report produced by this package is on an ExpressionSet
## of n arrays.  They are identified by the integers 1...n. I use the
## term 'report object' to be equivalent to 'array', and report
## objects are also identified by the same set of integers 1...n.
## This is in order to prepare the infrastructure created here also to
## be useable for reports on other types of objects than arrays (say:
## Solexa lanes, microtitre plates, ...).  The different plots in the
## report show points (e.g. in the PCA plot) or lines (e.g. in the
## density and so-called 'RNA-degradation' plots). Here, I call them
## 'plot objects'. They can be identified by the numbers
## 1...numPlotObjects .  In the simplest case, there is a 1:1
## correspondence between report objects and plot objects: there is
## exactly one line or point per array, and
## numPlotObjects==n. Sometimes however, we want to be more
## flexible. For instance, the density plot for two-color arrays shows
## *three* lines per array (one for red, one for green, one for the
## log-ratio), and numPlotObjects==3*n. Similarly, sometimes maybe we
## want to omit some arrays from a plot. To represent such mappings in
## general, the below functions are used.

setClass("svgParameters",
  representation(
    numPlotObjects      = "integer",             
    ## A function that finds the nodes in the SVG document corresponding to the plot objects.
    ## It should hold that length(getPlotObjNodes(doc)) == numPlotObjects            
    getPlotObjNodes     = "function",
    ## A function that given a report object ID (i.e. a number 1...n) returns the plot
    ## object ID (in the simplest case, this is a vector of length 1, but it can also
    ## have length 0 (if the object is not represented on this plot) or >1 (if there
    ## are more than 1 lines, points etc. for this report object)
    getPlotObjIdFromReportObjId  = "function",
    ## The inverse function of getPlotObjIdFromReportObjId
    getReportObjIdFromPlotObjId  = "function",   
    strokewidth         = "numeric",    ## vector of length 2: stroke-width without and with highlighting
    strokeopacity       = "numeric"),   ## vector of length 2: stroke-opacity without and with highlighting
         
  prototype(
    numPlotObjects      = NA_integer_,             
    getPlotObjNodes     = function(doc) NULL,
    getPlotObjIdFromReportObjId  = function(i) NA_integer_,
    getReportObjIdFromPlotObjId  = function(i) NA_integer_,
    strokewidth         = rep(NA_real_, 2),
    strokeopacity       = rep(NA_real_, 2)),
         
  validity = function(x) {
    if(length(x@numPlotObjects)!=1) return("Invalid slot 'numPlotObjects'.")
    if(length(x@strokewidth)   !=2) return("Invalid slot 'strokewidth'.")
    if(length(x@strokeopacity) !=2) return("Invalid slot 'strokeopacity'.")
    return(TRUE)
  }         
)


setClass("aqmReportModule",
  representation(
    plot     = "ANY",
    section  = "character",
    title    = "character",
    legend   = "character",
    shape    = "list",
    outliers = "integer",
    svg      = "list"),
  prototype(
    plot      = new("namedList"),
    section   = NA_character_,
    title     = NA_character_,
    legend    = NA_character_,
    shape     = list(),
    outliers  = integer(0),       
    svg       = list()))

setClassUnion("oneColourArray", c("ExpressionSet", "AffyBatch"))
