\name{primerIDAlignSeqs}
\alias{primerIDAlignSeqs}
\title{Align a short pattern with PrimerID to variable length target sequences.}
\usage{
primerIDAlignSeqs(subjectSeqs = NULL, patternSeq = NULL,
  qualityThreshold1 = 0.75, qualityThreshold2 = 0.5, doAnchored = FALSE,
  doRC = TRUE, returnUnmatched = FALSE, returnRejected = FALSE,
  showStats = FALSE, ...)
}
\arguments{
  \item{subjectSeqs}{DNAStringSet object containing
  sequences to be searched for the pattern.}

  \item{patternSeq}{DNAString object or a sequence
  containing the query sequence to search with the
  primerID.}

  \item{qualityThreshold1}{percent of first part of
  patternSeq to match. Default is 0.75.}

  \item{qualityThreshold2}{percent of second part of
  patternSeq to match. Default is 0.50.}

  \item{doAnchored}{for primerID based patternSeq, use the
  base before and after primer ID in patternSeq as
  anchors?. Default is FALSE.}

  \item{doRC}{perform reverse complement search of the
  defined pattern. Default is TRUE.}

  \item{returnUnmatched}{return sequences if it had no or
  less than 5\% match to the first part of patternSeq
  before the primerID. Default is FALSE.}

  \item{returnRejected}{return sequences if it only has a
  match to one side of patternSeq or primerID length does
  not match # of Ns +/-2 in the pattern. Default is FALSE.}

  \item{showStats}{toggle output of search statistics.
  Default is FALSE.}

  \item{...}{extra parameters for
  \code{\link{pairwiseAlignment}}}
}
\value{
\itemize{ \item A CompressedIRangesList of length two,
where x[["hits"]] is hits covering the entire patternSeq,
and x[["primerIDs"]] is the potential primerID region.
\item If returnUnmatched = T, then x[["Absent"]] is
appended which includes reads not matching the first part
of patternSeq.  \item If returnRejected=TRUE, then
x[["Rejected"]] includes reads that only matched one part
of patternSeq or places where no primerID was found in
between two part of patternSeq, and
x[["RejectedprimerIDs"]] includes primerIDs that didn't
match the correct length.  \item If doAnchored=TRUE, then
x[["unAnchoredprimerIDs"]] includes reads that didn't match
the base before and after primer ID on patternSeq. }
}
\description{
Align a fixed length short pattern sequence containing
primerID to variable length subject sequences using
\code{\link{pairwiseAlignment}}. This function uses default
of type="overlap", gapOpening=-1, and gapExtension=-1 to
align the patterSeq against subjectSeqs. The search is
broken up into as many pieces +1 as there are primerID and
then compared against subjectSeqs. For example,
patternSeq="AGCATCAGCANNNNNNNNNACGATCTACGCC" will launch
two search jobs one per either side of Ns. For each search,
qualityThreshold is used to filter out candidate alignments
and the area in between is chosen to be the primerID. This
strategy is benefical because of Indels introduced through
homopolymer errors. Most likely the length of primerID(s)
wont the same as you expected!
}
\note{
\itemize{ \item For qualityThreshold1 & qualityThreshold2,
the alignment score is calculated by
(matches*2)-(mismatches+gaps) which programatically
translates to round(nchar(patternSeq)*qualityThreshold)*2
\item Gaps and mismatches are weighed equally with value of
-1 which can be overriden by defining extra parameters
'gapOpening' & 'gapExtension'.  \item If qualityThreshold
is 1, then it is a full match, if 0, then any match is
accepted which is useful in searching linker sequences at
3' end. Beware, this function only searches for the pattern
sequence in one orientation. If you are expecting to find
the pattern in both orientation, you might be better off
using BLAST/BLAT! }
}
\examples{
#primerIDAlignSeqs(subjectSeqs, patternSeq, showStats=TRUE)
#primerIDAlignSeqs(subjectSeqs, patternSeq, showStats=TRUE, qualityThreshold1=0.5)
}
\seealso{
\code{\link{vpairwiseAlignSeqs}},
\code{\link{pairwiseAlignSeqs}}, \code{\link{doRCtest}},
\code{\link{blatSeqs}}, \code{\link{findAndRemoveVector}}
}

