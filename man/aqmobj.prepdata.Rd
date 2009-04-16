\name{aqmobj.prepdata-class}
\docType{class}
\alias{aqmobj.prepdata-class}
\alias{aqmobj.prepdata}
\alias{class.aqmobj.prepdata}

\title{ Class to contain data generated from aqm.prepdata. }

\description{Container for the output of \code{\link{aqm.prepdata}} and
  for the input of the aqm functions.}

\section{Slots}{
  \describe{
    \item{\code{M}:}{ A matrix of the M values (log-ratio). The
  log-ratio is computed with the second channel being the median of the
  intensities across arrays in the case of one channel arrays.}

    \item{\code{A}:}{ A matrix of the A values. The A value is the mean
  of the two intensities. The second channel is computed as for the M
  values in the case of one channel arrays.}

    \item{\code{dat}:}{ A matrix with the log-ratio if two channels or
      the intensities if one channel.}
    
    \item{\code{rc}:}{ A matrix with the red channel intensities in the
      case of two channels arrays.}
    
    \item{\code{gc}:}{ A matrix with the green channel intensities in the
      case of two channels arrays.}
 
    \item{\code{rcb}:}{ A matrix with the red channel background
      intensities if two channels arrays and if available.}
    
    \item{\code{gcb}:}{ A matrix with the green channel background
      intensities if two channels arrays and if available.}
    
    \item{\code{outM}:}{ The distance between each pairs of arrays,
      computed using \code{\link[genefilter:dist2]{dist2}} from the
      \code{\link[genefilter:dist2]{genefilter}} package.}
       
    \item{\code{sN}:}{ Integers numbering the arrays to be used to label
      the plots.}
    
    \item{\code{numArrays}:}{ An integer giving the number of arrays.}
    
    \item{\code{nchannels}:}{ A numeric giving the number of channels.}
    
    \item{\code{logtransformed}:}{ A logical telling if the data have
      been log transformed by the function aqm.prepdata.}
    
    \item{\code{classori}:}{ A character string of the class of the
      object that was given as an input of the aqm.prepdata function.}
    }
}
  
\details{
See the aqm.prepdata help or the aqm Vignette for example of this object.
}

\author{Audrey Kauffmann <audrey@ebi.ac.uk>}

\seealso{\code{\link{aqm.prepdata}}, \code{\link{aqm.boxplot}}, \code{\link{aqm.density}}, \code{\link{aqm.heatmap}}, \code{\link{aqm.maplot}}, \code{\link{aqm.meansd}}, \code{\link{aqm.probesmap}}, \code{\link{aqm.spatial}}, \code{\link{aqm.spatialbg}}}

\keyword{classes}
