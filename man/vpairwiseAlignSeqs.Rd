\name{vpairwiseAlignSeqs}
\alias{vpairwiseAlignSeqs}
\title{Align a pattern to variable length target sequences.}
\usage{
  vpairwiseAlignSeqs(subjectSeqs = NULL, patternSeq = NULL,
    side = "left", qualityThreshold = 1, showStats = FALSE,
    bufferBases = 5, doRC = TRUE, ...)
}
\arguments{
  \item{subjectSeqs}{DNAStringSet object containing
  sequences to be searched for the pattern. This is
  generally bigger than patternSeq, and cases where
  subjectSeqs is smaller than patternSeq will be ignored in
  the alignment.}

  \item{patternSeq}{DNAString object or a sequence
  containing the query sequence to search. This is
  generally smaller than subjectSeqs.}

  \item{side}{which side of the sequence to perform the
  search: left, right, or middle. Default is 'left'.}

  \item{qualityThreshold}{percent of patternSeq to match.
  Default is 1, full match. This is supplied to
  max.mismatch parameter of \code{\link{vmatchPattern}} as
  round(nchar(patternSeq)*(1-qualityThreshold)).}

  \item{showStats}{toggle output of search statistics.
  Default is FALSE.}

  \item{bufferBases}{use x number of bases in addition to
  patternSeq length to perform the search. Beneficial in
  cases where the pattern has homopolymers or indels
  compared to the subject. Default is 5. Doesn't apply when
  side='middle'.}

  \item{doRC}{perform reverse complement search of the
  defined pattern. Default is TRUE.}

  \item{...}{extra parameters for
  \code{\link{vmatchPattern}}}
}
\value{
  IRanges object with starts, stops, and names of the
  aligned sequences.
}
\description{
  Align a fixed length pattern sequence to subject
  sequences using \code{\link{vmatchPattern}}.
}
\note{
  Beware, this function only searches for the pattern
  sequence in one orientation. If you are expecting to find
  the pattern in both orientation, you might be better off
  using BLAST!
}
\examples{
#vpairwiseAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE)
 #vpairwiseAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE,qualityThreshold=0.5)
}
\seealso{
  \code{\link{pairwiseAlignSeqs}},
  \code{\link{primerIDAlignSeqs}}, \code{\link{doRCtest}},
  \code{\link{findAndTrimSeq}}
}

