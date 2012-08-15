\name{extractFeature}
\alias{extractFeature}
\title{Extract a specific feature/attribute of the sampleInfo object.
Given a sampleInfo object, the function extracts a defined feature(s) for given sample or sector.}
\usage{
  extractFeature(sampleInfo, sector = NULL,
    samplename = NULL, feature = NULL)
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object,
  which samples per sector/quadrant information along with
  other metadata.}

  \item{sector}{a vector or specific sector to extract the
  feature from. Default is NULL, which extracts all
  sectors.}

  \item{samplename}{a character vector or a specific sample
  to extract feature from. Default is NULL, which extracts
  all samples.}

  \item{feature}{Options include: primed, LTRed, linkered,
  decoded, and any of the metadata. Default is NULL. When
  feature='metadata', then it returns names of all the
  metadata elements associated with the sample as a comma
  separated list.}
}
\value{
  a list or list of lists depending upon which parameters
  were supplied.
}
\description{
  Extract a specific feature/attribute of the sampleInfo
  object. Given a sampleInfo object, the function extracts
  a defined feature(s) for given sample or sector.
}
\examples{
\dontrun{
 extractFeature(sampleInfo,feature="primed")
}
}
\seealso{
  \code{\link{addFeature}}, \code{\link{findPrimers}},
  \code{\link{findLTRs}}, \code{\link{findLinkers}},
  \code{\link{extractSeqs}}, \code{\link{trimSeqs}},
  \code{\link{getSectorsForSamples}}
}

