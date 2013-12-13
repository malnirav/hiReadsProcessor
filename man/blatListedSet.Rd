\name{blatListedSet}
\alias{blatListedSet}
\title{Align a listed DNAStringSet to a reference using gfClient or standalone BLAT.}
\usage{
blatListedSet(dnaSetList = NULL, ...)
}
\arguments{
  \item{dnaSetList}{DNAStringSet object containing
  sequences to be aligned against the reference.}

  \item{...}{parameters to be passed to
  \code{\link{blatSeqs}}.}
}
\value{
a list of RangedData object reflecting psl file type per
set of sequences.
}
\description{
Align sequences from a listed DNAStringSet object returned
from \code{\link{extractSeqs}} to an indexed reference
genome using gfServer/gfClient protocol or using standalone
BLAT and get a psl file as RangedData object. This function
heavily relies on defaults of \code{\link{blatSeqs}}.
}
\seealso{
\code{\link{pairwiseAlignSeqs}},
\code{\link{vpairwiseAlignSeqs}},
\code{\link{startgfServer}}, \code{\link{stopgfServer}},
\code{\link{blatSeqs}}, \code{\link{read.psl}},
\code{\link{pslToRangedObject}}, \code{\link{read.blast8}}
}

