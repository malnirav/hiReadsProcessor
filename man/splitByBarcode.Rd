\name{splitByBarcode}
\alias{splitByBarcode}
\title{Split DNAStringSet object using first X number of bases defined by a vector.
Given barcodes/MID to sample association vector and a DNAStringSet object, the function splits/demultiplexes the DNAStringSet object by first few bases dictated by length of barcodes/MID supplied. This is an accessory function used by \code{\link{decodeByBarcode}}}
\usage{
  splitByBarcode(barcodesSample, dnaSet, trimFrom = NULL,
    showStats = FALSE, returnUnmatched = FALSE,
    dereplicate = FALSE)
}
\arguments{
  \item{barcodesSample}{a character vector of barcodes to
  sample name associations. Ex: c("ACATCCAT"="Sample1",
  "GAATGGAT"="Sample2",...)}

  \item{dnaSet}{DNAStringSet object to evaluate.}

  \item{trimFrom}{integer value serving as start point to
  trim the sequences from. This is calculated internally
  length barcode+1. Default is NULL.}

  \item{showStats}{boolean flag denoting whether to show
  decoding statistics per sample & barcode. Default is
  FALSE.}

  \item{returnUnmatched}{boolean flag denoting whether to
  return unmatched reads. Default is FALSE.}

  \item{dereplicate}{return dereplicated sequences. Calls
  \code{\link{dereplicateReads}}, which appends counts=X to
  sequence names/deflines. Default is FALSE.}
}
\value{
  DNAStringSet object split by sample name found in
  barcodesSample.
}
\description{
  Split DNAStringSet object using first X number of bases
  defined by a vector. Given barcodes/MID to sample
  association vector and a DNAStringSet object, the
  function splits/demultiplexes the DNAStringSet object by
  first few bases dictated by length of barcodes/MID
  supplied. This is an accessory function used by
  \code{\link{decodeByBarcode}}
}
\examples{
\dontrun{
 splitByBarcode(c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2"),dnaSet,showStats=TRUE,returnUnmatched=TRUE)
}
}
\seealso{
  \code{\link{decodeByBarcode}},
  \code{\link{dereplicateReads}},
  \code{\link{replicateReads}}
}

