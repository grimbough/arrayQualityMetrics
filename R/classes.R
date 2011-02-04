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
## general, the below functions 'getPlotObjIdFromReportObjId' and
## 'getReportObjIdFromPlotObjId' are used.

setClass("svgParameters",
  representation(
    ## 'name' is used to construct the ID of Figures and Tables             
    name            = "character",
                 
    ## An R function that finds the nodes in the SVG document corresponding to the plot objects.
    ## It should hold that length(getPlotObjNodes(doc)) == numPlotObjects            
    getPlotObjNodes = "function",
    numPlotObjects  = "integer",
                 
    ## Two functions: The first, in JavaScript, when given a report
    ## object ID 'r' (e.g. "r:1") returns all associated plot object
    ## IDs 'p' (e.g. "p:1"). In the simplest case (which is what is
    ## this in the prototype definition below) this is a vector of
    ## length 1. However, it can also have length 0 (if the object is
    ## not represented on this plot) or >1 (if there are more than 1
    ## lines, points etc. for this report object). See the
    ## 'aqm.density' function for an example.
    ##
    ## The second function, in R, is the inverse. It is expected to
    ## always return a vector of length 1.
    getPlotObjIdFromReportObjId = "character",
    getReportObjIdFromPlotObjId = "function", 

    ## 2 x n matrix that contains stroke attributes such as stroke-width,
    ## stroke-opacity. Colnames of the matrix identify the attribute (e.g.
    ## colname 'opacity' corresponds to stroke-opacity). The matrix' first
    ## row corresponds to the value used for non-highlighted, the second
    ## row to those for highlighted objects.
    stroke     = "matrix"),   
         
  prototype(
    name            = NA_character_,
    getPlotObjNodes = getMatplotSeries,
    numPlotObjects  = NA_integer_,             
    getPlotObjIdFromReportObjId = "function(r) { return([r.replace(/^r:/, 'p:')]); }",
    getReportObjIdFromPlotObjId = function(x) sub("^p:", "r:", x),
    stroke          = matrix(c("1", "3", "0.4", "1"), nrow=2, dimnames = list(NULL, c("width", "opacity")))),
         
  validity = function(object) {
    if(length(object@name)          !=1) return("Invalid slot 'name'.")
    if(length(object@numPlotObjects)!=1) return("Invalid slot 'numPlotObjects'.")
    if(nrow(object@stroke)          !=2) return("Invalid slot 'stroke'.")
    return(TRUE)
  }         
)


##
## An object of this class contains everything needed to render a report module
##

setClass("aqmReportModule",
  representation(
    plot     = "ANY",
    size     = "numeric",     ## size of the plot in inch
    section  = "character",
    title    = "character",
    legend   = "character",
    outliers = "integer",
    svg      = "svgParameters"),

  prototype(
    plot      = new("namedList"),
    size      = c(w=NA_real_, h=NA_real_),
    section   = NA_character_,
    title     = NA_character_,
    legend    = NA_character_,
    outliers  = NA_integer_,       
    svg       = new("svgParameters")),

  validity = function(object) {
    for(s in c("section", "title", "legend"))
      if (length(slot(object, s)) != 1)
        return(sprintf("Invalid slot '%s'.", s))
    if ((length(object@size)!=2) || !identical(names(object@size), c("w", "h")))
      return("Invalid slot 'size'.")
    validObject(object@svg, test=TRUE)
  }
)

