\name{findIntegrations}
\alias{findIntegrations}
\title{Find the integration sites and add results to SampleInfo object.}
\usage{
  findIntegrations(sampleInfo, seqType = NULL, port = 5560,
    host = "localhost",
    genomeIndices = c(hg18 = "/usr/local/blatSuite34/hg18.2bit", mm8 = "/usr/local/blatSuite34/mm8.2bit"),
    parallel = TRUE, samplenames = NULL)
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object
  outputted from \code{\link{findLinkers}}, which holds
  decoded, primed, LTRed, and Linkered sequences for
  samples per sector/quadrant along with metadata.}

  \item{seqType}{which type of sequence to BLAT and find
  integration sites for. Default is NULL and determined
  automatically based on type of restriction enzyme or
  isolation method used. Could be any one of following
  options: genomic, genomicLinkered, decoded, primed,
  LTRed, linkered.}

  \item{port}{a port number to host the gfServer with.
  Default is 5560.}

  \item{host}{name of the machine running gfServer. Default
  is 'localhost'.}

  \item{genomeIndices}{an associative character vector of
  freeze to full or relative path of respective .nib or
  .2bit files. Default is
  c("hg18"="/usr/local/blatSuite34/hg18.2bit",
  "mm8"="/usr/local/blatSuite34/mm8.2bit").}

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
  sample. New data attributes include: psl, and sites. The
  psl attributes holds the BLAT hits per read along with QC
  information. The sites attribute holds the condensed
  integration sites where BLAT hits have been clustered by
  the Position column and cherry picked to have each site
  pass all the QC steps.
}
\description{
  Given a SampleInfo object, the function finds integration
  sites for each sample using their respective settings and
  adds the results back to the object. This is an
  all-in-one function which BLATs, parses PSL files, finds
  best hit per read per sample, cluster sites, and assign
  OTU IDs. It calls \code{\link{blatSeqs}},
  \code{\link{read.psl}},
  \code{\link{getIntegrationSites}},
  \code{\link{clusterSites}}, \code{\link{otuSites}}. There
  must be linkered reads within the sampleInfo object in
  order to use this function using the default parameters.
  If you are planning on BLATing non-linkered reads, then
  change the seqType to one of accepted options for the
  'feature' parameter of \code{\link{extractSeqs}}, except
  for '!' based features.
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
#findIntegrations(sampleInfo)
}
\seealso{
  \code{\link{findPrimers}}, \code{\link{findLTRs}},
  \code{\link{findLinkers}}, \code{\link{startgfServer}},
  \code{\link{read.psl}}, \code{\link{blatSeqs}},
  \code{\link{blatListedSet}},
  \code{\link{pslToRangedData}},
  \code{\link{clusterSites}}, \code{\link{otuSites}},
  \code{\link{getIntegrationSites}},
  \code{\link{read.blast8}}
}

