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
                    intgroup = NULL,
                    grouprep,
		    spatial = TRUE,
                    reporttitle = paste("Quality metrics report for", deparse(substitute(expressionset))))
}

\arguments{
  \item{expressionset}{an object of
    class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
    \code{\link[affy:AffyBatch-class]{AffyBatch}},
    \code{\link[Biobase:class.NChannelSet]{NChannelSet}},
    \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}},
    \code{\link[limma:rglist]{RGList}},
    \code{\link[limma:malist]{MAList}},
    \code{\link[marray:marrayRaw-class]{marrayRaw}},
    \code{\link[marray:marrayNorm-class]{marrayNorm}} or
    \code{aqmInputObj}}

\item{outdir}{the name of the directory in which the report is
  created; a character of length 1.}
  
\item{force}{if the directory named by \code{outdir} already exists, then, if \code{force} is
  \code{TRUE}, the directory is overwritten, otherwise an
  error is thrown; if the directory does not exist, the value of
  \code{force} is irrelevant; a logical of length 1.}
  
\item{do.logtransform}{indicates whether the data should be logarithm transformed before
  the analysis; a logical of length 1.}
  
\item{intgroup}{the name of the sample covariate(s) used to draw
  a colour side bar next to the heatmap. The first element of
  \code{intgroup} is also used define sample groups in
  other plots (boxplots, densities).
  \code{intgroup} should be a character vector, and its elements need to match the columns
  names of \code{phenoData(expressionset)}. If \code{NULL} or of length 0,
  then the plots are not decorated with sample covariate information.}

\item{grouprep}{deprecated. Use argument \code{intgroup} instead.}

\item{spatial}{indicates whether spatial plots
  should be made; a logical of length 1. This can be useful for large arrays (like Affymetrix
  hgu133Plus2) when CPU time and RAM resources of the machine would be limiting. }

\item{reporttitle}{title for the report (character of length 1).}
}

\details{
  See the arrayQualityMetrics vignette for examples of this function.
}

\value{
  A side effect of the function is the creation of directory named by \code{outdir} containing a HTML report \code{QMreport.html}
  and PNG and PDF files with figures. The function also returns a list with R
  objects containing the report elements for subsequent programmatic processing.
}

\author{
  Audrey Kauffmann.
}

\keyword{hplot}
