\name{read.psl}
\alias{read.psl}
\title{Read psl file(s) outputted by BLAT}
\usage{
  read.psl(pslFile = NULL, bestScoring = TRUE,
    asRangedData = FALSE, asGRanges = FALSE,
    removeFile = TRUE, parallel = FALSE)
}
\arguments{
  \item{pslFile}{psl filename, or vector of filenames, or a
  pattern of files to import.}

  \item{bestScoring}{report only best scoring hits instead
  of all hits. Default is TRUE. Score is calculated by
  matches-misMatches-qBaseInsert-tBaseInsert.}

  \item{asRangedData}{return a RangedData object instead of
  a dataframe. Default is FALSE}

  \item{asGRanges}{return a GRanges object instead of a
  dataframe. Default is FALSE}

  \item{removeFile}{remove the psl file(s) after importing.
  Default is FALSE.}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to TRUE.
  If no parallel backend is registered, then a serial
  version of foreach is ran using
  \code{\link{registerDoSEQ()}}.}
}
\value{
  a dataframe or RangedData object reflecting psl file
  type.
}
\description{
  Given filename(s), the function reads the psl file format
  from BLAT as a data frame and performs basic score
  filtering if indicated. Any other file format will yield
  errors or erroneous results.
}
\note{
  If parallel=TRUE, then be sure to have a paralle backend
  registered before running the function. One can use any
  of the following libraries compatible with
  \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For
  example: library(doMC); registerDoMC(2)
}
\examples{
#read.psl(pslFile="processed.*.psl$")
#read.psl(pslFile=c("sample1hits.psl","sample2hits.psl"))
}
\seealso{
  \code{\link{pairwiseAlignSeqs}},
  \code{\link{vpairwiseAlignSeqs}},
  \code{\link{startgfServer}}, \code{\link{blatSeqs}},
  \code{\link{read.blast8}}
}

