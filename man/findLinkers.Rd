\name{findLinkers}
\alias{findLinkers}
\title{Find the 3' linkers and add results to SampleInfo object.
Given a sampleInfo object, the function finds 3' linkers for each sample per sector and adds the results back to the object. This is a specialized function which depends on many other functions shown in 'see also section' to perform specialized trimming of 3' primer/linker adaptor sequence found in the sampleInfo object. The sequence itself is never trimmed but rather coordinates of linker portion is recorded back to the object and used subsequently by \code{\link{extractSeqs}} function to perform the trimming.}
\usage{
  findLinkers(sampleInfo, showStats = FALSE, doRC = FALSE,
    parallel = TRUE, samplenames = NULL)
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object
  outputted from \code{\link{findPrimers}} or
  \code{\link{findLTRs}}, which holds decoded sequences for
  samples per sector/quadrant along with information of
  sample to primer associations.}

  \item{showStats}{toggle output of search statistics.
  Default is FALSE.}

  \item{doRC}{perform reverse complement search of the
  defined pattern/linker sequence. Default is FALSE.}

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
  sample. New data attributes include: linkered. If linkers
  have primerID then, primerIDs attribute is appended as
  well.
}
\description{
  Find the 3' linkers and add results to SampleInfo object.
  Given a sampleInfo object, the function finds 3' linkers
  for each sample per sector and adds the results back to
  the object. This is a specialized function which depends
  on many other functions shown in 'see also section' to
  perform specialized trimming of 3' primer/linker adaptor
  sequence found in the sampleInfo object. The sequence
  itself is never trimmed but rather coordinates of linker
  portion is recorded back to the object and used
  subsequently by \code{\link{extractSeqs}} function to
  perform the trimming.
}
\note{
  If no linker matches are found with default options, then
  try doRC=TRUE. If parallel=TRUE, then be sure to have a
  paralle backend registered before running the function.
  One can use any of the following libraries compatible
  with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI.
  For example: library(doSMP); w <- startWorkers(2);
  registerDoSMP(w)
}
\examples{
\dontrun{
 findLinkers(sampleInfo,showStats=TRUE)
}
}
\seealso{
  \code{\link{pairwiseAlignSeqs}},
  \code{\link{vpairwiseAlignSeqs}},
  \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}},
  \code{\link{findPrimers}}, \code{\link{extractFeature}},
  \code{\link{extractSeqs}}, \code{\link{findAndTrimSeq}}
}

