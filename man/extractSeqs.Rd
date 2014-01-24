\name{extractSeqs}
\alias{extractSeqs}
\title{Extract sequences for a feature in the sampleInfo object.}
\usage{
extractSeqs(sampleInfo, sector = NULL, samplename = NULL,
  feature = "genomic", trim = TRUE, minReadLength = 1,
  sideReturn = NULL, pairReturn = "both")
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object,
  which samples per sector/quadrant information along with
  other metadata.}

  \item{sector}{specific sector to extract sequences from.
  Default is NULL, which extracts all sectors.}

  \item{samplename}{specific sample to extract sequences
  from. Default is NULL, which extracts all samples.}

  \item{feature}{which part of sequence to extract (case
  sensitive). Options include: primed, !primed, LTRed,
  !LTRed, linkered, !linkered, primerIDs, genomic,
  genomicLinkered, decoded, and unDecoded. If a sample was
  primerIDed and processed by
  \code{\link{primerIDAlignSeqs}}, then all the rejected
  and unmatched attributes can be prepended to the feature.
  Example: Rejectedlinkered, RejectedprimerIDslinkered,
  Absentlinkered, or unAnchoredprimerIDslinkered. When
  feature is genomic, it includes sequences which are
  primed, LTRed, linkered, and !linkered. The
  genomicLinkered is same as genomic minus the !linkered.
  When feature is decoded, it includes everything that
  demultiplexed. The '!' in front of a feature extracts the
  inverse. One can only get unDecoded sequences if
  returnUnmatched was TRUE in
  \code{\link{decodeByBarcode}}.}

  \item{trim}{whether to trim the given feature from
  sequences or keep it. Default is TRUE. This option is
  ignored for feature with '!'.}

  \item{minReadLength}{threshold for minimum length of
  trimmed sequences to return.}

  \item{sideReturn}{if trim=TRUE, which side of the
  sequence to return: left, middle, or right. Defaults to
  NULL and determined automatically. Doesn't apply to
  features: decoded, genomic or genomicLinkered.}

  \item{pairReturn}{if the data is paired end, then from
  which pair to return the feature. Options are "pair1",
  "pair2", or defaults to "both". Ignored if data is single
  end.}
}
\value{
a listed DNAStringSet object structed by sector then
sample. Note: when feature='genomic' or 'genomicLinkered'
and when data is paired end, then "pair2" includes union of
reads from both pairs which found LTR.
}
\description{
Given a sampleInfo object, the function extracts sequences
for a defined feature.
}
\examples{
#extractSeqs(sampleInfo)
#extractSeqs(sampleInfo,feature="primed")
}
\seealso{
\code{\link{findPrimers}}, \code{\link{findLTRs}},
\code{\link{findLinkers}}, \code{\link{trimSeqs}},
\code{\link{extractFeature}},
\code{\link{getSectorsForSamples}}
}

