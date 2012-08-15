\name{dereplicateReads}
\alias{dereplicateReads}
\title{Removes duplicate sequences from DNAStringSet object.
Given a DNAStringSet object, the function dereplicates reads and adds counts=X to the definition line to indicate replication.}
\usage{
  dereplicateReads(dnaSet)
}
\arguments{
  \item{dnaSet}{DNAStringSet object to dereplicate.}
}
\value{
  DNAStringSet object with names describing frequency of
  repeat.
}
\description{
  Removes duplicate sequences from DNAStringSet object.
  Given a DNAStringSet object, the function dereplicates
  reads and adds counts=X to the definition line to
  indicate replication.
}
\examples{
\dontrun{
 dereplicateReads(dnaSet)
}
}
\seealso{
  \code{\link{replicateReads}},
  \code{\link{removeReadsWithNs}},
  \code{\link{decodeByBarcode}},
  \code{\link{splitByBarcode}}
}

