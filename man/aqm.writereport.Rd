\name{aqm.writereport}
\alias{aqm.writereport}

\title{Write a quality report}

\description{
  \code{aqm.writereport} produces a quality report (HTML document with
  figures) from a list of \code{\link{aqmReportModule}} objects.}

\usage{aqm.writereport(name, obj, outdir) }

\arguments{
  \item{name, outdir}{Report title and output directory - as in \code{\link{arrayQualityMetrics}}.}
  \item{obj}{A list of \code{\link{aqmReportModule}} objects.}
}

\value{A side effect of this function is the creation of the HTML report.}

\details{
See the \code{aqm} Vignette.
}

\author{Audrey Kauffmann, Wolfgang Huber}
