## Terminology: I use the term 'report object' equivalent to 'array':
## this is in order to prepare the infrastructure created here also to
## be useable for reports on other types of objects than arrays (say:
## Solexa lanes, microtitre plates, ...).
##
## Suppose the report produced by this package is on an ExpressionSet
## of N_r arrays or 'report objects'.  They are identified by the
## identifiers ro:1, ro:2,... ro:(N_r).  The different plots in the
## report show points (e.g. in the PCA plot) or lines (e.g. in the
## density and so-called 'RNA-degradation' plots). Here, I call them
## 'plot objects'. They are identified by N_p identifiers po:1, po:2,
## ..., po:(N_p).  In the simplest case, there is a 1:1 correspondence
## between report objects and plot objects: there is exactly one line
## or point per array, and N_r=N_p.  Sometimes however, we want to be
## more flexible. For instance, the density plot for two-color arrays
## shows *three* lines per array (one for red, one for green, one for
## the log-ratio), and N_p=3*N_r. Similarly, sometimes maybe we want
## to omit some arrays from a plot. To represent such mappings in
## general, the below functions are used.
##
## Needed interactivity:
## - table refers to report objects. Checkbox update events receive
##   a report object id, and need to know how to update the related
##   plot object ids.
## - plots contain plot objects. Click and mouse over events receive
##   a plot object id, and need to know how to update the related
##   report object id(s) [for the table] and plot object ids.

setClass("svgParameters",
  representation(
    ## Objects of this class can defined (in which case the following slots
    ## contain useful information) or undefined.
    defined             = "logical",             
    numPlotObjects      = "integer",             
    ## An R function that finds the nodes in the SVG document corresponding to the plot objects.
    ## It should hold that length(getPlotObjNodes(doc)) == numPlotObjects            
    getPlotObjNodes     = "function",
    ## An array of two JavaScript functions:
    ## The first, when given a report object ID (e.g. "ro:1") returns all associated plot
    ## object IDs. The seccond is the inverse. The functions take a single string argument.
    ## The return a string array. In the simplest case, this is a vector of length 1, but it can also
    ## have length 0 (if the object is not represented on this plot) or >1 (if there
    ## are more than 1 lines, points etc. for this report object).
    idFun               = "character",
    tableID             = "character",
    strokewidth         = "numeric",    ## vector of length 2: stroke-width without and with highlighting
    strokeopacity       = "numeric"),   ## vector of length 2: stroke-opacity without and with highlighting
         
  prototype(
    defined          = FALSE,
    numPlotObjects   = NA_integer_,             
    getPlotObjNodes  = getMatplotSeries,
    idFun            = "[function(x) { [x.replace('^r', 'p')] }, function(x) { [x.replace('^p', 'r')] }]",
    tableID          = NA_character_,
    strokewidth      = c(1, 3),
    strokeopacity    = c(0.4, 1)),
         
  validity = function(x) {
    if(length(x@defined)!=1) return("Invalid slot 'defined'.")
    if(x@defined){
      if(length(x@numPlotObjects)!=1) return("Invalid slot 'numPlotObjects'.")
      if(length(x@tableID)       !=1) return("Invalid slot 'tableID'.")
      if(length(x@strokewidth)   !=2) return("Invalid slot 'strokewidth'.")
      if(length(x@strokeopacity) !=2) return("Invalid slot 'strokeopacity'.")
    }
    return(TRUE)
  }         
)


##
## All data needed to render a report Module
## 

setClass("aqmReportModule",
  representation(
    plot     = "ANY",
    section  = "character",
    title    = "character",
    legend   = "character",
    shape    = "list",
    outliers = "integer",
    svg      = "svgParameters"),

  prototype(
    plot      = new("namedList"),
    section   = NA_character_,
    title     = NA_character_,
    legend    = NA_character_,
    shape     = list(),
    outliers  = integer(0),       
    svg       = new("svgParameters")),

  validity = function(x) {
    for(s in c("section", "title", "legend"))
      if(length(slot(x, s))!=1) return(sprintf("Invalid slot '%s'.", s))
  }
)


##
## This is used by method definitions for 'prepdata' and 'spatial':
##

setClassUnion("oneColourArray", c("ExpressionSet", "AffyBatch"))
