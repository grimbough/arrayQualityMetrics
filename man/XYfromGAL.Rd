\name{addXYfromGAL}
\alias{addXYfromGAL-methods}
\alias{addXYfromGAL}

\title{
  Computing the coordinates of the spots on a slide 
}
\description{
  From the coordinates of the blocks of a microarray slide and the Row
  and Column locations of the spots within the blocks,
  \code{addXYfromGAL} computes the X and Y coordinates of the spots of a
  slide.  
}

\usage{
addXYfromGAL(x, gal.file, nBlocks, skip, ...)
}

\arguments{
  \item{x}{ is an \code{AnnotatedDataFrame} representing the featureData
    of an object. }
  
  \item{gal.file}{ name of the file .gal that contains the coordinates
    of the blocks. }
  
  \item{nBlocks}{ number of blocks on the slide. }
  
  \item{skip}{ number of header lines to skip when reading the
    gal.file. }
  
  \item{...}{Arguments that get passed on to \code{read.table}.}  
}

\value{
  The object \code{x} of class \code{AnnotatedDataFrame} will be
  returned with two added columns: X and Y corresponding to the absolute
  position of the probes on the array.
  }

\author{
Audrey Kauffmann, Wolfgang Huber.
Maintainer: <kauffmann@bergonie.org>
}
\keyword{file}
\keyword{manip}
