\name{findAndRemoveVector}
\alias{findAndRemoveVector}
\title{Find and trim vector sequence from reads.}
\usage{
findAndRemoveVector(reads, Vector, minLength = 10, returnCoords = FALSE,
  parallel = TRUE)
}
\arguments{
  \item{reads}{DNAStringSet object containing sequences to
  be trimmed for vector.}

  \item{Vector}{DNAString object containing vector sequence
  to be searched in reads.}

  \item{minLength}{integer value dictating minimum length
  of trimmed sequences to return. Default is 10.}

  \item{returnCoords}{return the coordinates of vector
  start-stop for the matching reads. Defaults to FALSE.}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to TRUE.
  If no parallel backend is registered, then a serial
  version of foreach is ran using
  \code{\link{registerDoSEQ()}}.}
}
\value{
DNAStringSet object with Vector sequence removed from the
reads. If returnCoords=TRUE, then a list of two named
elements "hits" & "reads". The first element, "hits" is a
GRanges object with properties of matched region and
whether it was considered valid denoted by 'good.row'. The
second element, "reads" is a DNAStringSet object with
Vector sequence removed from the reads.
}
\description{
This function facilitates finding and trimming of
long/short fragments of vector present in LM-PCR products.
The algorithm looks for vector sequence present anywhere
within the read and trims according longest contiguous
match on either end of the read. Alignment is doing using
BLAT
}
\note{
If parallel=TRUE, then be sure to have a paralle backend
registered before running the function. One can use any of
the following libraries compatible with
\code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For
example: library(doMC); registerDoMC(2)
}
\examples{
#findAndRemoveVector(reads, Vector)
}
\seealso{
\code{\link{pairwiseAlignSeqs}},
\code{\link{vpairwiseAlignSeqs}},
\code{\link{pslToRangedObject}}, \code{\link{blatSeqs}},
\code{\link{read.blast8}}, \code{\link{findAndTrimSeq}}
}

