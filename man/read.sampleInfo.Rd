\name{read.sampleInfo}
\alias{read.sampleInfo}
\title{Read a sample information file and format appropriate metadata.}
\usage{
  read.sampleInfo(sampleInfoPath = NULL,
    splitBySector = TRUE, interactive = TRUE)
}
\arguments{
  \item{sampleInfoFilePath}{full or relative path to the
  sample information file, which holds samples to
  quadrant/lane associations along with other metadata
  required to trim sequences or process it.}

  \item{splitBySector}{split the data frame into a list by
  sector column. Default is TRUE.}

  \item{interactive}{whether to prompt each time the
  function encounters an issue, or use the defaults.
  Default is TRUE.}
}
\value{
  if splitBySector=TRUE, then an object of SimpleList named
  by quadrant/lane information defined in sampleInfo file,
  else a dataframe.
}
\description{
  Given a sample information file, the function checks if
  it includes required information to process samples
  present on each sector/quadrant/region/lane. The function
  also adds other columns required for processing with
  default values if not already defined ahead of time.
}
\details{
  \itemize{ \item Required Column Description: \itemize{
  \item sector => region/quadrant of the sequencing plate
  the sample comes from. If files have been split by
  samples apriori, then the filename associated per sample
  without the extension. If this is a filename, then be
  sure to enable 'alreadyDecoded' parameter in
  \code{\link{findBarcodes}} or
  \code{\link{decodeByBarcode}}, since contents of this
  column is pasted together with 'seqfilePattern' parameter
  in \code{\link{read.SeqFolder}} to find the appropriate
  file needed.  \item barcode => unique 4-12bp DNA sequence
  which identifies the sample. If providing filename as
  sector, then leave this blank.  \item primerltrsequence
  => DNA sequence of the viral LTR primer with viral LTR
  sequence following the primer landing site. If already
  trimmed, then mark this as SKIP.  \item samplename =>
  Name of the sample associated with the barcode \item
  sampledescription => Detailed description of the sample
  \item gender => sex of the sample: male or female \item
  species => species of the sample: homo sapien, mus
  musculus, etc.  \item freeze => UCSC freeze to which the
  sample should be aligned to.  \item linkersequence => DNA
  sequence of linker adaptor to found following the genomic
  sequence. If already trimmed, then mark this as SKIP.
  \item restrictionenzyme => Restriction enzyme used for
  digestion and sample recovery. Can also be Fragmentase or
  Sonication!  } \item Metadata Parameter Column
  Description: \itemize{ \item ltrbitsequence => DNA
  sequence of the viral LTR following the primer landing
  site.  \item ltrbitidentity => percent of LTR bit
  sequence to match during the alignment. Default is 1.
  \item primerltridentity => percent of primer to match
  during the alignment. Default is .85 \item linkeridentity
  => percent of linker sequence to match during the
  alignment. Default is 0.55. Only applies to
  non-primerID/random tag based linker search.  \item
  primeridinlinker => whether the linker adaptor used has
  primerID/random tag in it? Default is FALSE.  \item
  primeridinlinkeridentity1 => percent of sequence to match
  before the random tag. Default is 0.75. Only applies to
  primerID/random tag based linker search and when
  primeridinlinker is TRUE.  \item
  primeridinlinkeridentity2 => percent of sequence to match
  after the random tag. Default is 0.50. Only applies to
  primerID/random tag based linker search and when
  primeridinlinker is TRUE.  \item celltype => celltype
  information associated with the sample \item user => name
  of the user who prepared or processed the sample } \item
  Processing Parameter Column Description: \itemize{ \item
  startwithin => upper bound limit of where the alignment
  should start within the query. Default is 3.  \item
  alignratiothreshold => cuttoff for (alignment span/read
  length). Default is 0.7.  \item genomicpercentidentity =>
  cuttoff for (1-(misMatches/matches)). Default is 0.98.
  \item clustersiteswithin => cluster integration sites
  within a defined window size based on frequency which
  corrects for any sequencing errors. Default is 5.  \item
  keepmultihits => whether to keep sequences/reads that
  return multiple best hits, aka ambiguous locations.
  \item processingdate => the date of processing } }
}
\examples{
#read.sampleInfo("~/sampleInfo.xls")
}
\seealso{
  \code{\link{read.SeqFolder}},
  \code{\link{decodeByBarcode}},
  \code{\link{splitByBarcode}}
}

