\name{read.BAMasPSL}
\alias{read.BAMasPSL}
\title{Reads a BAM/SAM file and converts it into a PSL like format.}
\usage{
read.BAMasPSL(bamFile = NULL, removeFile = TRUE, asGRanges = TRUE)
}
\arguments{
  \item{bamFile}{BAM/SAM filename, or vector of filenames,
  or a pattern of files to import.}

  \item{removeFile}{remove the file(s) supplied in bamFile
  paramter after importing. Default is FALSE.}

  \item{asGRanges}{return a GRanges object. Default is
  TRUE}
}
\value{
a GRanges or GAlignments object reflecting psl file type.
}
\description{
Given filename(s), the function reads the BAM/SAM file,
converts into a PSL like format. Any other file format will
yield errors or erroneous results.
}
\examples{
#read.BAMasPSL(bamFile="processed.*.bam$")
#read.BAMasPSL(bamFile=c("sample1hits.sam","sample2hits.sam"))
}
\seealso{
\code{\link{pairwiseAlignSeqs}}, \code{\link{blatSeqs}},
\code{\link{read.blast8}}, \code{\link{read.psl}},
\code{\link{pslToRangedObject}},
\code{\link{subreadAlignSeqs}},
\code{\link{pairUpAlignments}}
}

