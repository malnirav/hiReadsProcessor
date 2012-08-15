\name{removeReadsWithNs}
\alias{removeReadsWithNs}
\title{Remove sequences with ambiguious nucleotides.
Given a DNAStringSet object, the function removes any reads that has either repeating or total Ns which is greater than to maxNs threshold}
\usage{
  removeReadsWithNs(dnaSet, maxNs = 5, consecutive = TRUE)
}
\arguments{
  \item{dnaSet}{DNAStringSet object to evaluate.}

  \item{maxNs}{integer value denoting the threshold of
  maximum allowed Ns. Default is 5.}

  \item{consecutive}{boolean flag denoting whether Ns to
  filter is consecutive or total . Default is TRUE.}
}
\value{
  DNAStringSet object.
}
\description{
  Remove sequences with ambiguious nucleotides. Given a
  DNAStringSet object, the function removes any reads that
  has either repeating or total Ns which is greater than to
  maxNs threshold
}
\examples{
\dontrun{
 removeReadsWithNs(dnaSet)
}
}
\seealso{
  \code{\link{dereplicateReads}},
  \code{\link{replicateReads}},
  \code{\link{decodeByBarcode}},
  \code{\link{splitByBarcode}}
}

