\name{arrayQualityMetrics}
\alias{arrayQualityMetrics}
\alias{arrayQualityMetrics-methods}
\alias{arrayQualityMetrics,ExpressionSet-method}
\alias{arrayQualityMetrics,AffyBatch-method}
\alias{arrayQualityMetrics,NChannelSet-method}

\title{
  Quality metrics on microarray experiments
}

\description{
  \code{arrayQualityMetrics} performs quality metrics on
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
  \code{\link[Biobase:class.AffyBatch]{AffyBatch}} or
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}} containing
  microarray data from any platforms, one or two channels. The results,
  presented in a HTML report, are designated to allow the user
  to rapidly assess the quality of a set of arrays.
}

\usage{
\S4method{arrayQualityMetrics}{ExpressionSet}(expressionset, outdir, force, do.logtransform, split.plots)

arrayQualityMetrics(expressionset,
                    outdir = getwd(),
                    force = FALSE,
                    do.logtransform = FALSE,
                    split.plots = FALSE)

\S4method{arrayQualityMetrics}{AffyBatch}(expressionset, outdir, force, do.logtransform, split.plots)
\S4method{arrayQualityMetrics}{NChannelSet}(expressionset, outdir, force, do.logtransform, split.plots)
}

\arguments{
  \item{expressionset}{ is an object of
    class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
    \code{\link[Biobase:class.AffyBatch]{AffyBatch}} or
    \code{\link[Biobase:class.NChannelSet]{NChannelSet}}.}

\item{outdir}{ is the name of the directory in which the results are
  created. }
  
\item{force}{ if TRUE, \code{outdir} will be overwritten if it
  already exists. }
  
\item{do.logtransform}{ If TRUE, the data are log transformed before
  the analysis. }
  
\item{split.plots}{  If the number of studied array
  is more than 50 it is adviced to define a number of experiments to
  represent on the density plots. }
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
                        do.logtransform = TRUE,
                        split.plots = 10)
    # "MLL_QMreport.html" is a report created in the working directory
    # on the logarithm expression data of the MLL.A AffyBatch.
    # 10 experiments are represented by density plots.
    }
}

\keyword{hplot}
