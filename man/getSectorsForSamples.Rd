\name{getSectorsForSamples}
\alias{getSectorsForSamples}
\title{Get sectors for samples defined in the sampleInfo object.}
\usage{
  getSectorsForSamples(sampleInfo, sector = NULL,
    samplename = NULL, returnDf = FALSE)
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object,
  which samples per sector/quadrant information along with
  other metadata.}

  \item{sector}{a specific sector or vector of sectors if
  known ahead of time. Default is NULL, which extracts all
  sectors.}

  \item{samplename}{a specific sample or vector of
  samplenames to get sectors for. Default is NULL, which
  extracts all samples.}

  \item{returnDf}{return results in a dataframe. Default is
  FALSE.}
}
\value{
  If returnDf=TRUE, then a dataframe of sector associated
  with each samplename, else a named list of length two:
  x[["sectors"]] and x[["samplenames"]]
}
\description{
  Given a sampleInfo object, the function gets the sectors
  for each samplename. This is an accessory function
  utilized by other functions of this package to aid sector
  retrieval.
}
\examples{
#getSectorsForSamples(sampleInfo,samplename="SampleName1")
}
\seealso{
  \code{\link{extractSeqs}}, \code{\link{extractFeature}},
  \code{\link{addFeature}}
}

