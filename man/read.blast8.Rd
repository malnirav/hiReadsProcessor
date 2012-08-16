\name{read.blast8}
\alias{read.blast8}
\title{Read blast8 file(s) outputted by BLAT}
\usage{
  read.blast8(files = NULL, asRangedData = FALSE,
    removeFile = TRUE, parallel = FALSE)
}
\arguments{
  \item{files}{blast8 filename, or vector of filenames, or
  a pattern of files to import.}

  \item{asRangedData}{return a RangedData object instead of
  a dataframe. Default is TRUE. Saves memory!}

  \item{removeFile}{remove the blast8 file(s) after
  importing. Default is FALSE.}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to TRUE.
  If no parallel backend is registered, then a serial
  version of foreach is ran using
  \code{\link{registerDoSEQ()}}.}
}
\value{
  a dataframe or RangedData object reflecting blast8 file
  type.
}
\description{
  Given filename(s), the function reads the blast8 file
  format from BLAT as a data frame and performs basic score
  filtering if indicated. Any other file format will yield
  errors or erroneous results.
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
#read.blast8(files="processed.*.blast8$")
 #read.blast8(files=c("sample1hits.blast8","sample2hits.blast8"))
}
\seealso{
  \code{\link{pairwiseAlignSeqs}},
  \code{\link{vpairwiseAlignSeqs}},
  \code{\link{startgfServer}}, \code{\link{blatSeqs}}
}

