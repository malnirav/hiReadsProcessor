\name{blatSeqs}
\alias{blatSeqs}
\title{Align sequences using BLAT.}
\usage{
  blatSeqs(query = NULL, subject = NULL,
    standaloneBlat = TRUE, port = NULL, host = "localhost",
    parallel = TRUE, gzipResults = TRUE,
    blatParameters = c(minIdentity = 70, minScore = 5, stepSize = 5, tileSize = 10, repMatch = 112312, dots = 50, q = "dna", t = "dna", out = "psl"))
}
\arguments{
  \item{query}{an object of DNAStringSet, a character
  vector, or a path/pattern of fasta files to BLAT. Default
  is NULL.}

  \item{subject}{an object of DNAStringSet, a character
  vector, or a path to an indexed genome (nibs,2bits) to
  serve as a reference or target to the query. Default is
  NULL. If the subject is a path to a nib or 2bit file,
  then standaloneBlat will not work!}

  \item{standaloneBlat}{use standalone BLAT as suppose to
  gfServer/gfClient protocol. Default is TRUE.}

  \item{port}{the same number you started the gfServer
  with. Required if standaloneBlat=FALSE. Default is NULL.}

  \item{host}{name of the machine running gfServer. Default
  is 'localhost' and only used when standaloneBlat=FALSE.}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to TRUE.
  If no parallel backend is registered, then a serial
  version of foreach is ran using
  \code{\link{registerDoSEQ()}}.}

  \item{gzipResults}{gzip the output files? Default is
  TRUE.}

  \item{blatParameters}{a character vector of options to be
  passed to gfClient/BLAT command except for 'nohead'
  option. Default: c(minIdentity=70, minScore=5,
  stepSize=5, tileSize=10, repMatch=112312, dots=50,
  q="dna", t="dna", out="psl"). Be sure to only pass
  parameters accepted by either BLAT or gfClient. For
  example, if repMatch or stepSize parameters are specified
  when using gfClient, then the function will simply ignore
  them! The defaults are configured to align a 19bp
  sequence with 70% identity.}
}
\value{
  a character vector of psl filenames.
}
\description{
  Align batch of sequences using standalone BLAT or
  gfServer/gfClient protocol for alignment against an
  indexed reference genome. Depending on parameters
  provided, the function either aligns batch of files to a
  reference genome using gfClient or takes sequences from
  query & subject parameters and aligns them using
  standalone BLAT. If standaloneBlat=FALSE and gfServer is
  not launched apriori, this function will start one using
  \code{\link{startgfServer}} and kill it using
  \code{\link{stopgfServer}} upon successful execution.
}
\examples{
\dontrun{
 blatSeqs(dnaSeqs,subjectSeqs,blatParameters=c(minIdentity=90, minScore=10, tileSize=10, dots=10, q="dna", t="dna", out="blast8"))
 blatSeqs(dnaSeqs,"/usr/local/blatSuite34/hg18.2bit",standaloneBlat=FALSE)
}
}
\seealso{
  \code{\link{pairwiseAlignSeqs}},
  \code{\link{vpairwiseAlignSeqs}},
  \code{\link{startgfServer}}, \code{\link{stopgfServer}},
  \code{\link{read.psl}}, \code{\link{splitSeqsToFiles}}
}

