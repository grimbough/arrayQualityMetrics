\name{aqm.writereport}
\alias{aqm.writereport}

\title{ Write HTML reports from objects produced with aqm functions. }

\description{
  \code{aqm.writereport} performs a HTML report including summary with outliers detected, titles, plots and legends.}

\usage{ aqm.writereport(expressionset, dataprep, intgroup = "Covariate", ...) }

\arguments{
  \item{name}{ A name to customize the title of the report that will be "name quality metrics report" }
  \item{expressionset}{ The expressionset on which the metrics have been ran }
  \item{obj}{ A list of aqmobj.xxx objects }
}

\value{ A HTML report "QMreport.html" in the working directory. }

\author{ Audrey Kauffmann <audrey@ebi.ac.uk> }

\examples{
\dontrun{
}
}