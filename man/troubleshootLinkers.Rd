\name{troubleshootLinkers}
\alias{troubleshootLinkers}
\title{Compare LTRed sequences to all linkers.}
\usage{
  troubleshootLinkers(sampleInfo, qualityThreshold = 0.55,
    qualityThreshold1 = 0.75, qualityThreshold2 = 0.5,
    doRC = TRUE, parallel = TRUE, samplenames = NULL)
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object
  outputted from \code{\link{findPrimers}} or
  \code{\link{findLTRs}}, which holds decoded sequences for
  samples per sector/quadrant along with information of
  sample to primer associations.}

  \item{qualityThreshold}{percent of linker length to
  match, round(nchar(linker)*qualityThreshold). Default is
  0.55. Only applies to non-primerID based linkers}

  \item{qualityThreshold1}{percent of first part of
  patternSeq to match. Default is 0.75. Only applies to
  primerID based linker search.}

  \item{qualityThreshold2}{percent of second part of
  patternSeq to match. Default is 0.50. Only applies to
  primerID based linker search.}

  \item{doRC}{perform reverse complement search of the
  linker sequence. Default is TRUE. Highly recommended!}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to TRUE.
  If no parallel backend is registered, then a serial
  version of foreach is ran using
  \code{\link{registerDoSEQ()}}.}

  \item{samplenames}{a vector of samplenames to process.
  Default is NULL, which processes all samples from
  sampleInfo object.}
}
\value{
  a dataframe of counts.
}
\description{
  Given a SampleInfo object, the function compares LTRed
  sequences from each sample per sector to all the linker
  sequences present in the run. The output is a summary
  table of counts of good matches to all the linkers per
  sample.
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
#troubleshootLinkers(sampleInfo,showStats=TRUE)
}
\seealso{
  \code{\link{pairwiseAlignSeqs}},
  \code{\link{vpairwiseAlignSeqs}},
  \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}},
  \code{\link{findPrimers}}, \code{\link{findAndTrimSeq}}
}

