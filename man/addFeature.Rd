\name{addFeature}
\alias{addFeature}
\title{Add a specific feature/attribute to the sampleInfo object.
Given a sampleInfo object, the function adds a new feature for the given samples & sectors.}
\usage{
  addFeature(sampleInfo, sector = NULL, samplename = NULL,
    feature = NULL, value = NULL)
}
\arguments{
  \item{sampleInfo}{sample information SimpleList object,
  which samples per sector/quadrant information along with
  other metadata.}

  \item{sector}{a vector or a specific sector to add the
  new feature(s) to. Default is NULL, in which case the
  sectors are searched from samplename parameter.}

  \item{samplename}{a character vector or a specific sample
  to add the new feature(s) to. Default is NULL.}

  \item{feature}{a string of naming the new feature to add
  for the defined samplename and sector.}

  \item{value}{the value or a named list of samplenames &
  values which is assigned for the defined sector,
  samplename, and feature. Example:
  list("Sample1"="ACDTDASD")}
}
\value{
  modified sampleInfo object with new feature(s) added.
}
\description{
  Add a specific feature/attribute to the sampleInfo
  object. Given a sampleInfo object, the function adds a
  new feature for the given samples & sectors.
}
\examples{
\dontrun{
 addFeature(sampleInfo,"1","Sample1",feature="primerIDed","ACDTDASD")
}
}
\seealso{
  \code{\link{findPrimers}}, \code{\link{extractSeqs}},
  \code{\link{trimSeqs}}, \code{\link{extractFeature}},
  \code{\link{getSectorsForSamples}}
}

