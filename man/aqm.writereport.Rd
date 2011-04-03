\name{aqm.writereport}
\alias{aqm.writereport}

\title{Write a quality report}

\description{
  \code{aqm.writereport} produces a quality report (HTML document with
  figures) from a list of \code{\link{aqmReportModule}} objects.}

\usage{aqm.writereport(modules, arrayTable, reporttitle, outdir)}

\arguments{
  \item{modules}{A list of \code{\link{aqmReportModule}} objects.}
  \item{arrayTable}{A data.frame with array (meta)data to be displayed
  in the report.}
  \item{reporttitle, outdir}{Report title and output directory - as in \code{\link{arrayQualityMetrics}}.}
}

\value{A side effect of this function is the creation of the HTML report.}

\author{Audrey Kauffmann, Wolfgang Huber}
