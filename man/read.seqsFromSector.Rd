\name{read.seqsFromSector}
\alias{read.seqsFromSector}
\title{Read fasta/fastq given the path or sampleInfo object.}
\usage{
read.seqsFromSector(seqFilePath = NULL, sector = 1)
}
\arguments{
  \item{seqFilePath}{a path to fasta/fastq reads file or a
  sampleInfo object returned by
  \code{\link{read.SeqFolder}}}

  \item{sector}{specific sector to reads sequences from.
  Default is 1, and not required if seqFilePath is a direct
  file path rather than sampleInfo object.}
}
\value{
a DNAStringSet object.
}
\description{
Given a sequence reads file path, the function returns a
DNAStringSet object.
}
\seealso{
\code{\link{decodeByBarcode}},
\code{\link{read.SeqFolder}}, \code{\link{extractSeqs}}
}

