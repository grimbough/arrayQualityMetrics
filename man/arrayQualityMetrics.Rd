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

\title{ Quality metrics on microarray experiments }

\description{
  \code{arrayQualityMetrics} performs quality metrics on
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
                    intgroup = "Covariate",
                    grouprep = FALSE,
		    spatial = TRUE)

}

\arguments{
  \item{expressionset}{ is an object of
    class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
    \code{\link[affy:AffyBatch-class]{AffyBatch}} ,
    \code{\link[Biobase:class.NChannelSet]{NChannelSet}} ,
    \code{\link[beadarray:class-BeadLevelList]{BeadLevelList}} ,
    \code{\link[limma:rglist]{RGList}} ,
    \code{\link[limma:malist]{MAList}} ,
    \code{aqmInputObj} ,
    \code{\link[marray:marrayRaw-class]{marrayRaw}} or
    \code{\link[marray:marrayNorm-class]{marrayNorm}} .}

\item{outdir}{ is the name of the directory in which the results are
  created. }
  
\item{force}{ if TRUE, \code{outdir} will be overwritten if it
  already exists. }
  
\item{do.logtransform}{ If TRUE, the data are log transformed before
  the analysis. }
  
\item{intgroup}{ Name of the column of the phenoData to be used to draw
  a colour side bar next to the heatmap. }

\item{grouprep}{ Decide if you want the boxplots and density plots to be
  coloured function of the groups set by 'intgroup'. The default
  is FALSE meaning that the boxplot and density plots will not be
  represented function of the groups of 'intgroup'.}

\item{spatial}{ If FALSE, the spatial representations are not
  performed. This is useful for large arrays(like Affymetrix
  hgu133Plus2) when these figures are slow to be drawn and can cause
  machines with low memory to fail to perform the report. } 
}

\details{
See the arrayQualityMetrics Vignette for examples of this function.
}

\value{
  A directory \code{outdir} containing a HTML report named QMreport.html
  and all the PNG and PDF plots is created. 
}

\author{
  Audrey Kauffmann, Wolfgang Huber. 
  Maintainer: <audrey@ebi.ac.uk>  
}

\keyword{hplot}
