\name{arrayQualityMetrics-package}
\alias{arrayQualityMetrics}
\docType{package}
\title{
Quality metrics on ExpressionSets
}
\description{
  This package performs quality metrics on AffyBatch, ExpressionSet
  or NChannelSet containing  microarray data from any platforms,
  one or two channels. The results are designated to allow the user
  to rapidly assess the quality of a set of arrays.
}
\usage{
arrayQualityMetrics(expressionset,
                    outdir = getwd(),
                    force = FALSE,
                    do.logtransform = FALSE,
                    split.plots = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expressionset}{ is an object of
    class ExpressionSet, AffyBatch or NChannelSet. }
  \item{outdir}{ is the name of the directory in which the results are
    created. }
    \item{force}{ if TRUE, if outdir already exists, it will be overwritten. }
  \item{do.logtransform}{ If TRUE, the data are log transformed before
    the analysis. }
  \item{split.plots}{  If the number of studied array
    is more than 50 it is adviced to define a number of experiments to
    represent on the density plots. }
}
\details{
\tabular{ll}{
Package: \tab arrayQualityMetrics\cr
Type: \tab Package\cr
Version: \tab 1.0.13\cr
Date: \tab 2007-08-06\cr
License: \tab LGPL \cr
}
}
\author{
Audrey Kauffmann, Wolfgang Huber,
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
\keyword{ hplot }
\keyword{ graphs }% __ONLY ONE__ keyword per line
