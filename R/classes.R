## Terminology: I use the term 'report object' equivalent to 'array'.
## This is in order to prepare the infrastructure created here also to
## be useable for reports on other types of objects than arrays (say:
## Solexa lanes, microtitre plates, ...). Report objects are identified 
## by the numbers 1, 2,... N_r.
##
setClass("svgParameters",
  representation(
    ## An R function that finds the nodes in the SVG document corresponding to the plot objects.
    ## It should hold that length(getPlotObjNodes(doc)) == numPlotObjects            
    getPlotObjNodes = "function",
    numPlotObjects  = "integer",
    getReportObjIdFromPlotObjId = "function"),   
         
  prototype(
    getPlotObjNodes = getMatplotSeries,
    numPlotObjects  = NA_integer_,             
    getReportObjIdFromPlotObjId = function(x) x),
         
  validity = function(object) {
    if(length(object@numPlotObjects)!=1) return("Invalid slot 'numPlotObjects'.")
    return(TRUE)
  }         
)

##
## Parameters summarising an outlier detection
##
setClass("outlierDetection",
  representation(
    statistic   = "numeric",
    threshold   = "numeric",
    which       = "integer",
    description = "character"),   
  prototype(
    statistic   = numeric(0),
    threshold   = NA_real_,
    which       = NA_integer_,
    description = NA_character_))

##
## An object of this class contains everything needed to render a report module
##
setClass("aqmReportModule",
  representation(
    plot           = "ANY",
    size           = "numeric",     ## size of the plot in inch
    colors         = "character",
    section        = "character",
    title          = "character",
    id             = "character",
    legend         = "character",
    outliers       = "outlierDetection",
    defaultdisplay = "character",             
    svg            = "svgParameters"),

  prototype(
    plot           = new("namedList"),
    size           = c(w=NA_real_, h=NA_real_),
    colors         = "#b0b0b0",
    section        = NA_character_,
    title          = NA_character_,
    id             = NA_character_,
    legend         = NA_character_,
    outliers       = new("outlierDetection"),
    defaultdisplay = "block",      
    svg            = new("svgParameters")),

  validity = function(object) {
    for(s in c("section", "title", "legend"))
      if (length(slot(object, s)) != 1)
        return(sprintf("Invalid slot '%s'.", s))
    if ((length(object@size)!=2) || !identical(names(object@size), c("w", "h")))
      return("Invalid slot 'size'.")
    if(!(is.character(object@defaultdisplay) && (length(object@defaultdisplay)==1) && (object@defaultdisplay %in% c("block", "none"))))
      return("Invalid slot 'defaultdisplay'.")
    validObject(object@svg, test=TRUE)
  }
)

