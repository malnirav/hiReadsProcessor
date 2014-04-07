\name{pairUpAlignments}
\alias{pairUpAlignments}
\title{Pair up alignments in a GRanges object}
\usage{
pairUpAlignments(psl.rd = NULL, maxGapLength = 2500, sameStrand = TRUE,
  parallel = TRUE)
}
\arguments{
  \item{psl.rd}{a GRanges object with qNames ending in "}

  \item{maxGapLength}{maximum gap allowed between end of
  pair1 and start of pair2. Default is 1000 bp.}

  \item{sameStrand}{should pairs be aligned to the same
  strand or in same orientationin the reference genome?
  Default is TRUE. This is 'TRUE' because pair2 reads are
  reverseComplemented when reading in data in
  \code{\link{findBarcodes}}}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to TRUE.
  If no parallel backend is registered, then a serial
  version of foreach is ran using
  \code{\link{registerDoSEQ()}}.}
}
\value{
a GRanges object with reads paired up denoted by "paired"
column. Improper pairs or unpaired reads are returned with
"paired" column as FALSE.
}
\description{
Given a GRanges object, the function uses specified
gaplength parameter to pair up reads where the qName column
ends with "
}
\examples{
#pairUpAlignments()
#pairUpAlignments()
}
\seealso{
\code{\link{pairwiseAlignSeqs}}, \code{\link{blatSeqs}},
\code{\link{read.blast8}}, \code{\link{read.psl}},
\code{\link{getIntegrationSites}},
\code{\link{subreadAlignSeqs}}, \code{\link{read.BAMasPSL}}
}

