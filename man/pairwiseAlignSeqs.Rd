\name{pairwiseAlignSeqs}
\alias{pairwiseAlignSeqs}
\title{Align a pattern to variable length target sequences.}
\usage{
  pairwiseAlignSeqs(subjectSeqs = NULL, patternSeq = NULL,
    side = "left", qualityThreshold = 1, showStats = FALSE,
    bufferBases = 5, doRC = TRUE, returnUnmatched = FALSE,
    returnLowScored = FALSE, ...)
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
  search: left, right or middle. Default is 'left'.}

  \item{qualityThreshold}{percent of patternSeq to match.
  Default is 1, full match.}

  \item{showStats}{toggle output of search statistics.
  Default is FALSE.}

  \item{bufferBases}{use x number of bases in addition to
  patternSeq length to perform the search. Beneficial in
  cases where the pattern has homopolymers or indels
  compared to the subject. Default is 5. Doesn't apply when
  side='middle'.}

  \item{doRC}{perform reverse complement search of the
  defined pattern. Default is TRUE.}

  \item{returnUnmatched}{return sequences which had no or
  less than 5\% match to the patternSeq. Default is FALSE.}

  \item{returnLowScored}{return sequences which had quality
  score less than the defined qualityThreshold. Default is
  FALSE.}

  \item{...}{extra parameters for
  \code{\link{pairwiseAlignment}}}
}
\value{
  IRanges object with starts, stops, and names of the
  aligned sequences. If returnLowScored or returnUnmatched
  = T, then a CompressedIRangesList where x[["hits"]] has
  the good scoring hits, x[["Rejected"]] has the failed to
  match qualityThreshold hits, and x[["Absent"]] has the
  hits where the aligned bit is <=10\% match to the
  patternSeq.
}
\description{
  Align a fixed length pattern sequence to subject
  sequences using \code{\link{pairwiseAlignment}}. This
  function uses default of type="overlap", gapOpening=-1,
  and gapExtension=-1 to align the patternSeq against
  subjectSeqs. One can adjust these parameters if prefered,
  but not recommended.
}
\note{
  For qualityThreshold, the alignment score is calculated
  by (matches*2)-(mismatches+gaps) which programatically
  translates to
  round(nchar(patternSeq)*qualityThreshold)*2. Gaps and
  mismatches are weighed equally with value of -1. If
  qualityThreshold is 1, then its is a full match, if 0,
  then any match is accepted which is useful in searching
  linker sequences at 3' end. Beware, this function only
  searches for the pattern sequence in one orientation. If
  you are expecting to find the pattern in both
  orientation, you might be better off using BLAST!
}
\examples{
#pairwiseAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE)
 #pairwiseAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE,qualityThreshold=0.5)
}
\seealso{
  \code{\link{primerIDAlignSeqs}},
  \code{\link{vpairwiseAlignSeqs}}, \code{\link{doRCtest}},
  \code{\link{findAndTrimSeq}}
}

