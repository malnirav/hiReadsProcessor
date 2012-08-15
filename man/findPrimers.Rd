\name{findPrimers}
\alias{findPrimers}
\title{Find the 5' primers and add results to SampleInfo object.}
\usage{
  findPrimers(sampleInfo, alignWay = "slow",
    showStats = FALSE, doRC = FALSE, parallel = TRUE,
    samplenames = NULL)
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object
  outputted from \code{\link{decodeByBarcode}}, which holds
  decoded sequences for samples per sector/quadrant along
  with information of sample to primer associations.}

  \item{alignWay}{method to utilize for detecting the
  primers. One of following: "slow" (Default), or "fast".
  Fast, calls \code{\link{vpairwiseAlignSeqs}} and uses
  \code{\link{vpatternMatch}} at its core, which is less
  accurate with indels and mismatches but much faster.
  Slow, calls \code{\link{pairwiseAlignSeqs}} and uses
  \code{\link{pairwiseAlignment}} at its core, which is
  accurate with indels and mismatches but slower.}

  \item{showStats}{toggle output of search statistics.
  Default is FALSE.}

  \item{doRC}{perform reverse complement search of the
  defined pattern/primer. Default is FALSE.}

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
  a SimpleList object similar to sampleInfo paramter
  supplied with new data added under each sector and
  sample. New data attributes include: primed
}
\description{
  Given a sampleInfo object, the function finds 5' primers
  for each sample per sector and adds the results back to
  the object. This is a specialized function which depends
  on many other functions shown in 'see also section' to
  perform specialized trimming of 5' primer/adaptor found
  in the sampleInfo object. The sequence itself is never
  trimmed but rather coordinates of primer portion is
  recorded back to the object and used subsequently by
  \code{\link{extractSeqs}} function to perform the
  trimming.
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
 findPrimers(sampleInfo,showStats=TRUE)
 findPrimers(sampleInfo,alignWay="slow",showStats=TRUE)
}
}
\seealso{
  \code{\link{pairwiseAlignSeqs}},
  \code{\link{vpairwiseAlignSeqs}},
  \code{\link{extractFeature}}, \code{\link{extractSeqs}},
  \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}},
  \code{\link{findLinkers}}, \code{\link{findAndTrimSeq}}
}

