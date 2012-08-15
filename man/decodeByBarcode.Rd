\name{decodeByBarcode}
\alias{decodeByBarcode}
\title{Demultiplex reads by their barcodes}
\usage{
  decodeByBarcode(sampleInfo, sector = NULL, dnaSet = NULL,
    showStats = FALSE, returnUnmatched = FALSE,
    dereplicate = FALSE)
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object
  created using \code{\link{read.sampleInfo}}, which holds
  barcodes and sample names per sector/quadrant or a
  character vector of barcodes to sample name associations.
  Ex: c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2",...)}

  \item{sector}{If sampleInfo is a SimpleList object, then
  a numeric/character value or vector representing
  sector(s) in sampleInfo. Optionally if on high memory
  machine sector='all' will decode/demultiplex sequences
  from all sectors/quadrants. This option is ignored if
  sampleInfo is a character vector. Default is NULL.}

  \item{dnaSet}{DNAStringSet object containing sequences to
  be decoded or demultiplexed. Default is NULL. If
  sampleInfo is a SimpleList object, then reads are
  automatically extracted using
  \code{\link{read.seqsFromSector}} and parameters defined
  in sampleInfo object.}

  \item{showStats}{toggle output of search statistics.
  Default is FALSE.}

  \item{returnUnmatched}{return unmatched sequences.
  Returns results as a list where x[["unDecodedSeqs"]] has
  culprits. Default is FALSE.}

  \item{dereplicate}{return dereplicated sequences. Calls
  \code{\link{dereplicateReads}}, which appends counts=X to
  sequence names/deflines. Default is FALSE.}
}
\value{
  If sampleInfo is an object of SimpleList then decoded
  sequences are appeneded to respective sample slots, else
  a named list of DNAStringSet object. If
  returnUnmatched=TRUE, then x[["unDecodedSeqs"]] has the
  unmatched sequences.
}
\description{
  Given a sample information object, the function reads in
  the raw fasta/fastq file, demultiplexes reads by their
  barcodes, and appends it back to the sampleInfo object.
  Calls \code{\link{splitByBarcode}} to perform the actual
  splitting. If supplied with a character vector and reads
  themselves, the function behaves a bit differently. See
  the examples.
}
\examples{
\dontrun{
 decodeByBarcode(sampleInfo,showStats=TRUE)
 decodeByBarcode(sampleInfo=c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2"),dnaSet,showStats=TRUE,returnUnmatched=TRUE)
}
}
\seealso{
  \code{\link{splitByBarcode}},
  \code{\link{dereplicateReads}},
  \code{\link{replicateReads}}
}

