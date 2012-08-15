\name{splitSeqsToFiles}
\alias{splitSeqsToFiles}
\title{Split DNA sequences into smaller files.}
\usage{
  splitSeqsToFiles(x, totalFiles = 4, suffix = "tempy",
    filename = "queryFile.fa")
}
\arguments{
  \item{x}{a DNAStringSet object, or a FASTA filename.}

  \item{totalFiles}{an integer indicating how many files to
  create. Default is 4.}

  \item{suffix}{a word to add to each file created. Default
  is "tempy".}

  \item{filename}{name of the file if x is a DNAStringSet
  object. Default is "queryFile.fa".}
}
\value{
  a vector of filename names created.
}
\description{
  Given a vector of sequences or DNAStringSet or a FASTA
  filename, the function splits it into smaller pieces as
  denoted by totalFiles parameter.
}
\examples{
\dontrun{
 splitSeqsToFiles(dnaSeq,10,"tempyQ","myDNAseqs.fa")
 splitSeqsToFiles("mySeqs.fa",5,"tempyQ")
}
}
\seealso{
  \code{\link{blatSeqs}}
}

