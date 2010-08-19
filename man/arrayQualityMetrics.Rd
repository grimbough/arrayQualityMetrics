\name{arrayQualityMetrics}
\alias{arrayQualityMetrics}
\alias{arrayQualityMetrics-methods}
\alias{arrayQualityMetrics,ExpressionSet-method}
\alias{arrayQualityMetrics,AffyBatch-method}
\alias{arrayQualityMetrics,NChannelSet-method}
\alias{arrayQualityMetrics,BeadLevelList-method}
\alias{arrayQualityMetrics,RGList-method}
\alias{arrayQualityMetrics,MAList-method}
\alias{arrayQualityMetrics,marrayRaw-method}
\alias{arrayQualityMetrics,marrayNorm-method}
\alias{arrayQualityMetrics,aqmInputObj-method}

\title{Quality metrics for microarray experiments}

\description{
  \code{arrayQualityMetrics} reports quality metrics on
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
  \code{\link[affy:AffyBatch-class]{AffyBatch}} ,
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}},
  \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}},
  \code{\link[limma:rglist]{RGList}},
  \code{\link[limma:malist]{MAList}},
  \code{aqmInputObj},
  \code{\link[marray:marrayRaw-class]{marrayRaw}} or
  \code{\link[marray:marrayNorm-class]{marrayNorm}} containing
  microarray data from any platforms, one or two channels. The results,
  presented in a HTML report, are designated to allow the user
  to rapidly assess the quality of a set of arrays.
}

\usage{
arrayQualityMetrics(expressionset,
                    outdir = getwd(),
                    force = FALSE,
                    do.logtransform = FALSE,
                    intgroup,
                    grouprep,
		    spatial = TRUE,
		    sN = NULL,
                    reporttitle = paste("Quality metrics report for", deparse(substitute(expressionset))))
}

\arguments{
  \item{expressionset}{an object of
    class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
    \code{\link[affy:AffyBatch-class]{AffyBatch}} ,
    \code{\link[Biobase:class.NChannelSet]{NChannelSet}} ,
    \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}} ,
    \code{\link[limma:rglist]{RGList}} ,
    \code{\link[limma:malist]{MAList}} ,
    \code{aqmInputObj} ,
    \code{\link[marray:marrayRaw-class]{marrayRaw}} or
    \code{\link[marray:marrayNorm-class]{marrayNorm}} .}

\item{outdir}{the name of the directory in which the results are
  created; a character of length 1.}
  
\item{force}{indicates whether \code{outdir} should be overwritten if it
  already exists; a logical of length 1.}
  
\item{do.logtransform}{indicates whether the data should be logarithm transformed before
  the analysis; a logical of length 1.}
  
\item{intgroup}{the name of the sample covariate(s) used to draw
  a colour side bar next to the heatmap. The first element of
  \code{intgroup} is also used define sample groups in
  other plots (boxplots, densities).
  \code{intgroup} should be a character vector, and its elements need to match the columns
  names of \code{phenoData(expressionset)}. If missing or \code{NA},
  then the plots are not decorated with sample covariate information.}

\item{grouprep}{Deprecated. This argument used to indicate whether or not plots should
  be decorated with sample covariate information (see argument \code{intgroup}).}

\item{spatial}{indicates whether spatial plots
  should be made; a logical of length 1. This can be useful for large arrays (like Affymetrix
  hgu133Plus2) when CPU time and RAM resources of the machine would be limiting. }

\item{sN}{Should the sample names be written on the plots? If
  \code{NULL}, samples will be labeled by numbers 1, 2, .... on the
  plots, and a table is produced to map these numbers back to the
  names.}
\item{reporttitle}{Title for the report (character of length 1).}
}

\details{
See the arrayQualityMetrics vignette for examples of this function.
}

\value{
  A directory \code{outdir} containing a HTML report named QMreport.html
  and all the PNG and PDF plots is created. 
}

\author{
  Audrey Kauffmann, Wolfgang Huber. 
  Maintainer: <kauffmann@bergonie.org>  
}

\keyword{hplot}
