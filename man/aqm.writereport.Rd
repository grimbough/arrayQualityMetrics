\name{aqm.writereport}
\alias{aqm.writereport}

\title{ Writes a report from objects produced with aqm.xxx functions. }

\description{
  \code{aqm.writereport} performs an html report from a list
  of aqmobj objects. It includes a summary with the outliers detected, titles, plots and legends.}

\usage{ aqm.writereport(name, expressionset, obj) }

\arguments{
  \item{name}{ A name to customize the title of the report that will be "name quality metrics report" }
  \item{expressionset}{ The expressionset on which the metrics have been run }
  \item{obj}{ A list of aqmobj.xxx objects }
}

\value{ An html report named 'QMreport.html' in the working directory. }

\details{
See the aqm Vignette for example of this object.
}

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }
