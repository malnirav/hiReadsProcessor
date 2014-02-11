\name{findIntegrations}
\alias{findIntegrations}
\title{Find the integration sites and add results to SampleInfo object.}
\usage{
findIntegrations(sampleInfo, seqType = NULL, aligner = "BLAT",
  port = 5560, host = "localhost", genomeIndices = NULL,
  parallel = TRUE, samplenames = NULL, blatParameters = c(minIdentity =
  70, minScore = 5, stepSize = 5, tileSize = 10, repMatch = 112312, dots = 50, q
  = "dna", t = "dna", out = "psl"))
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object
  outputted from \code{\link{findLinkers}}, which holds
  decoded, primed, LTRed, and Linkered sequences for
  samples per sector/quadrant along with metadata.}

  \item{seqType}{which type of sequence to align and find
  integration sites for. Default is NULL and determined
  automatically based on type of restriction enzyme or
  isolation method used. If restriction enzyme is
  Fragmentase or Mu, then this parameter is set to
  genomicLinkered, else it is genomic. Any one of following
  options are valid: genomic, genomicLinkered, decoded,
  primed, LTRed, linkered.}

  \item{aligner}{which aligner to use: BLAT(default) or
  Subread. Use 'Subread' for illumina data.}

  \item{port}{a port number to host the gfServer with. Only
  used if aligner='BLAT'. Default is 5560.}

  \item{host}{name of the machine running gfServer. Only
  used if aligner='BLAT'. Default is 'localhost'.}

  \item{genomeIndices}{an associative character vector of
  freeze to full or relative path of respective of indexed
  genomes from BLAT(.nib or .2bit files) or
  Subreads(.subread.index). For example:
  c("hg18"="/usr/local/blatSuite34/hg18.2bit",
  "mm8"="/usr/local/blatSuite34/mm8.2bit"). Be sure to
  supply an index per freeze supplied in the sampleInfo
  object. Default is NULL.}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to TRUE.
  If no parallel backend is registered, then a serial
  version of foreach is ran using
  \code{\link{registerDoSEQ()}}.}

  \item{samplenames}{a vector of samplenames to process.
  Default is NULL, which processes all samples from
  sampleInfo object.}

  \item{blatParameters}{a character vector of options to be
  passed to gfClient/BLAT command except for 'nohead'
  option. Default: c(minIdentity=70, minScore=5,
  stepSize=5, tileSize=10, repMatch=112312, dots=50,
  q="dna", t="dna", out="psl"). Be sure to only pass
  parameters accepted by either BLAT or gfClient. For
  example, if repMatch or stepSize parameters are specified
  when using gfClient, then the function will simply ignore
  them! The defaults are configured to align a 19bp
  sequence with 70\% identity. Ignored if aligner is
  "Subread".}
}
\value{
a SimpleList object similar to sampleInfo parameter
supplied with new data added under each sector and sample.
New data attributes include: psl, and sites. The psl
attributes holds the genomic hits per read along with QC
information. The sites attribute holds the condensed
integration sites where genomic hits have been clustered by
the Position column and cherry picked to have each site
pass all the QC steps.
}
\description{
Given a SampleInfo object, the function finds integration
sites for each sample using their respective settings and
adds the results back to the object. This is an all-in-one
function which aligns, finds best hit per read per sample,
cluster sites, and assign ISU IDs. Depending on the aligner
chosen, it calls \code{\link{blatSeqs}},
\code{\link{read.psl}}, \code{\link{getIntegrationSites}},
\code{\link{clusterSites}}, \code{\link{otuSites}}. There
must be linkered reads within the sampleInfo object in
order to use this function using the default parameters. If
you are planning on BLATing non-linkered reads, then change
the seqType to one of accepted options for the 'feature'
parameter of \code{\link{extractSeqs}}, except for '!'
based features.
}
\note{
If parallel=TRUE, then be sure to have a paralle backend
registered before running the function. One can use any of
the following libraries compatible with
\code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For
example: library(doMC); registerDoMC(2)
}
\examples{
#findIntegrations(sampleInfo)
}
\seealso{
\code{\link{findPrimers}}, \code{\link{findLTRs}},
\code{\link{findLinkers}}, \code{\link{startgfServer}},
\code{\link{read.psl}}, \code{\link{blatSeqs}},
\code{\link{blatListedSet}},
\code{\link{pslToRangedObject}},
\code{\link{clusterSites}}, \code{\link{otuSites2}},
\code{\link{crossOverCheck}},
\code{\link{getIntegrationSites}}
}

