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
  \code{\link[affy:class.AffyBatch]{AffyBatch}} ,
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}},
  \code{\link[beadarray:class.BeadLevelList]{BeadLevelList}},
  \code{\link[limma:class.RGList]{RGList}},
  \code{\link[limma:class.MAList]{MAList}},
  \code{\link[limma:class.aqmInputObj]{aqmInputObj}},
  \code{\link[marray:class.marrayRaw]{marrayRaw}} or
  \code{\link[marray:class.marrayNorm]{marrayNorm}} containing
  microarray data from any platforms, one or two channels. The results,
  presented in a HTML report, are designated to allow the user
  to rapidly assess the quality of a set of arrays.
}

\usage{
\S4method{arrayQualityMetrics}{ExpressionSet}(expressionset, outdir,
force, do.logtransform, split.plots, intgroup, grouprep)

arrayQualityMetrics(expressionset,
                    outdir = getwd(),
                    force = FALSE,
                    do.logtransform = FALSE,
                    split.plots = FALSE,
                    intgroup = "Covariate",
                    grouprep = FALSE)

\S4method{arrayQualityMetrics}{AffyBatch}(expressionset, outdir, force,
do.logtransform, split.plots, intgroup, grouprep)
\S4method{arrayQualityMetrics}{NChannelSet}(expressionset, outdir,
force, do.logtransform, split.plots, intgroup, grouprep)
\S4method{arrayQualityMetrics}{BeadLevelList}(expressionset, outdir,
force, do.logtransform, split.plots, intgroup, grouprep)
\S4method{arrayQualityMetrics}{RGList}(expressionset, outdir,
force, do.logtransform, split.plots, intgroup, grouprep)
\S4method{arrayQualityMetrics}{MAList}(expressionset, outdir,
force, do.logtransform, split.plots, intgroup, grouprep)
\S4method{arrayQualityMetrics}{aqmInputObj}(expressionset, outdir,
force, do.logtransform, split.plots, intgroup, grouprep)
\S4method{arrayQualityMetrics}{marrayRaw}(expressionset, outdir,
force, do.logtransform, split.plots, intgroup, grouprep)
\S4method{arrayQualityMetrics}{marrayNorm}(expressionset, outdir,
force, do.logtransform, split.plots, intgroup, grouprep)
}

\arguments{
  \item{expressionset}{ is an object of
    class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
    \code{\link[Biobase:class.AffyBatch]{AffyBatch}} ,
    \code{\link[Biobase:class.NChannelSet]{NChannelSet}} ,
    \code{\link[Biobase:class.BeadLevelList]{BeadLevelList}} ,
    \code{\link[Biobase:class.RGList]{RGList}} ,
    \code{\link[Biobase:class.MAList]{MAList}} ,
    \code{\link[Biobase:class.aqmInputObj]{aqmInputObj}} ,
    \code{\link[Biobase:class.marrayRaw]{marrayRaw}} or
    \code{\link[Biobase:class.marrayNorm]{marrayNorm}} .}

\item{outdir}{ is the name of the directory in which the results are
  created. }
  
\item{force}{ if TRUE, \code{outdir} will be overwritten if it
  already exists. }
  
\item{do.logtransform}{ If TRUE, the data are log transformed before
  the analysis. }
  
\item{split.plots}{ If the number of studied array
  is more than 50 it is adviced to define a number of experiments to
  represent on the density plots. }

\item{intgroup}{ Name of the column of the phenoData to be used to draw
  a color side bar next to the heatmap. }

\item{grouprep}{ Decide if you want the boxplots and density plots to be
  coloured function of the groups set by 'intgroup'. The default
  is FALSE meaning that the boxplot and density plots will not be
  represented function of the groups of 'intgroup'.}
}

\value{
  A directory \code{outdir} containing a HTML report named QMreport.html and all the PNG and PDF plots is created.
}

\author{
  Audrey Kauffmann, Wolfgang Huber. 
  Maintainer: <audrey@ebi.ac.uk>  
}

\examples{
  \dontrun{
    library(ALLMLL)
    data(MLL.A)
    arrayQualityMetrics(expressionset = MLL.A,
                        outdir ="MLL",
                        do.logtransform = TRUE)
    # "QMreport.html" is a report created in the "MLL" directory
    # on the logarithm expression data of the MLL.A AffyBatch.
    }
}

\keyword{hplot}
