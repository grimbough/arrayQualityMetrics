\name{aqm.writereport}
\alias{aqm.writereport}

\title{ Writes a report from objects produced with aqm.xxx functions. }

\description{
  \code{aqm.writereport} performs an html report from a list
  of aqmobj objects. It includes a summary with the outliers detected, titles, plots and legends.}

\usage{aqm.writereport(name, obj, outdir) }

\arguments{
  \item{name}{Report title.}
  \item{obj}{A list of objects with report sections.}
  \item{obj}{ A list of aqmobj.xxx objects }
}

\value{A side effect of this function is the creation of the HTML report.}

\details{
See the \code{aqm} Vignette.
}

\author{Audrey Kauffmann}
