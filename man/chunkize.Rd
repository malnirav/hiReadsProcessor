\name{chunkize}
\alias{chunkize}
\title{Breaks an object into chunks of N size.}
\usage{
chunkize(x, chunkSize = NULL)
}
\arguments{
  \item{x}{a linear object.}

  \item{chunkSize}{number of rows to use per chunk of
  query. Default to length(sites.rd)/detectCores() or
  length(query)/getDoParWorkers() depending on parallel
  backend registered.}
}
\value{
a list of object split into chunks.
}
\description{
Given a linear object, the function breaks query into
chunks of N size where each chunk has a respective subject
object filtered by seqnames/space present in the query
chunk. This is a helper function used by functions in 'See
Also' section where each chunk is sent to a parallel node
for processing.
}
\examples{
chunkize(DNAStringSet(c("GAGGCTGTCACCGACAAGGTTCTGA", "AATAGCGTGGTGACAGCCCACATGC", "GGTCTTCTAGGGAACCTACGCCACA",
"TTTCCGGCGGCAGTCAGAGCCAAAG", "TCCTGTCAACTCGTAGATCCAATCA", "ATGGTCACCTACACACAACGGCAAT",
"GTCAGGACACCTAATCACAAACGGA", "AGACGCAGGTTCAGGCTGGACTAGA", "ATCGTTTCCGGAATTCGTGCACTGG",
"CAATGCGGGCACACGCTCTTACAGT")), 5)
}
\seealso{
\code{\link{primerIDAlignSeqs}},
\code{\link{vpairwiseAlignSeqs}},
\code{\link{pairwiseAlignSeqs}}
}

