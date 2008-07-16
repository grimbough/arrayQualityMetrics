\name{rg2ncs}
\alias{rg2ncs}

\title{
  Convert a RGList into a NChannelSet}
\description{
  From an object of class RGList, \code{rg2ncs} builds a minimal
  NChannelSet with the expression values and the featureNames.
}

\usage{
rg2ncs(RG, type = "complete")
}

\arguments{
  \item{RG}{ is an \code{RGList}. }
  
  \item{type}{ "complete" if the RGList contains the slots R, G, Rb and Gb.
    "nobkg" if the RGList contains the slots R and G.
  }
}
\value{
  An NChannelSet with assayData and featureNames. }

\author{
Audrey Kauffmann.
Maintainer: <audrey@ebi.ac.uk>
}
\keyword{file}
\keyword{manip}
