\name{findAndTrimSeq}
\alias{findAndTrimSeq}
\title{Find and trim the pattern sequence from the subject.}
\usage{
  findAndTrimSeq(patternSeq, subjectSeqs, side = "left",
    offBy = 0, alignWay = "slow", ...)
}
\arguments{
  \item{patternSeq}{DNAString object or a sequence
  containing the query sequence to search.}

  \item{subjectSeqs}{DNAStringSet object containing
  sequences to be searched for the pattern.}

  \item{side}{which side of the sequence to perform the
  search & trimming: left, right or middle. Default is
  'left'.}

  \item{offBy}{integer value dictating if the trimming base
  should be offset by X number of bases. Default is 0.}

  \item{alignWay}{method to utilize for detecting the
  primers. One of following: "slow" (Default), or "fast".
  Fast, calls \code{\link{vpairwiseAlignSeqs}} and uses
  \code{\link{vpatternMatch}} at its core, which is less
  accurate with indels and mismatches but much faster.
  Slow, calls \code{\link{pairwiseAlignSeqs}} and uses
  \code{\link{pairwiseAlignment}} at its core, which is
  accurate with indels and mismatches but slower.}

  \item{...}{parameters to be passed to
  \code{\link{pairwiseAlignment}} or
  \code{\link{vpairwiseAlignSeqs}} depending on which
  method is defined in 'alignWay' parameter.}
}
\value{
  DNAStringSet object with pattern sequence removed from
  the subject sequences.
}
\description{
  This function facilitates finding and trimming of the
  pattern sequence from a collection of subject sequences.
  The trimming is dictated by side parameter. For more
  information on the trimming process see the 'side'
  parameter documentation in \code{\link{trimSeqs}}. For
  information regarding the pattern alignment see the
  documentation for \code{\link{pairwiseAlignSeqs}}.
}
\note{
  If parallel=TRUE, then be sure to have a paralle backend
  registered before running the function. One can use any
  of the following libraries compatible with
  \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For
  example: library(doSMP); w <- startWorkers(2);
  registerDoSMP(w)
}
\examples{
\dontrun{
findAndTrimSeq(patternSeq="AGACCCTTTT",subjectSeqs=DNAStringSet(c("AGACCCTTTTGAGCAGCAT","AGACCCTTGGTCGACTCA","AGACCCTTTTGACGAGCTAG")), qualityThreshold=.85, doRC=F, side="left", offBy=1, alignWay = "slow")
}
}
\seealso{
  \code{\link{pairwiseAlignSeqs}},
  \code{\link{vpairwiseAlignSeqs}},
  \code{\link{extractFeature}}, \code{\link{extractSeqs}},
  \code{\link{primerIDAlignSeqs}},
  \code{\link{findPrimers}}, \code{\link{findLinkers}}
}

