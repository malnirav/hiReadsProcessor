\name{otuSites2}
\alias{otuSites2}
\title{Bin values or make OTUs by assigning a unique ID to them within discrete factors.}
\usage{
  otuSites2(posID = NULL, value = NULL, readID = NULL,
    grouping = NULL, psl.rd = NULL, parallel = TRUE)
}
\arguments{
  \item{posID}{a vector of groupings for the value
  parameter (i.e. Chr,strand). Required if psl.rd parameter
  is not defined.}

  \item{value}{a vector of integer locations/positions that
  needs to be binned, i.e. genomic location. Required if
  psl.rd parameter is not defined.}

  \item{readID}{a vector of read/clone names which is
  unique to each row, i.e. deflines.}

  \item{grouping}{additional vector of grouping by which to
  pool the rows (i.e. samplenames). Default is NULL.}

  \item{psl.rd}{a RangedData/GRanges object returned from
  \code{\link{clusterSites}}. Default is NULL.}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to TRUE.
  If no parallel backend is registered, then a serial
  version of foreach is ran using
  \code{\link{registerDoSEQ()}}. Process is split by the
  grouping the column.}
}
\value{
  a data frame with binned values and otuID shown alongside
  the original input. If psl.rd parameter is defined, then
  a RangedData/GRanges object where object is first
  filtered by clusterTopHit column and the otuID column
  appended at the end.
}
\description{
  Given a group of values or genomic positions per
  read/clone, the function tries to yield a unique OTU
  (operation taxinomical unit) ID for the collection based
  on overlap of locations to other reads/clones by
  grouping. This is mainly useful when each read has many
  locations which needs to be considered as one single
  group of sites.
}
\note{
  The algorithm for making OTUs of sites is as follows: for
  each readID check how many positions are there. Separate
  readIDs with only position from the rest. Check if any
  readIDs with >1 position match to any readIDs with only
  one position. If there is a match, then assign both
  readIDs with the same OTU ID. Check if any positions from
  readIDs with >1 position match any other readIDs with >1
  position. If yes, then assign same OTU ID to all readIDs
  sharing 1 or more positions.
}
\examples{
#otuSites2(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), value=c(rep(1000,2),5832,1000,12324,65738,928042), readID=paste('read',sample(letters,7),sep='-'), grouping=c('a','a','a','b','b','b','c'))
#otuSites2(psl.rd=test.psl.rd)
}
\seealso{
  \code{\link{clusterSites}}, \code{\link{otuSites}},
  \code{\link{crossOverCheck}},
  \code{\link{findIntegrations}},
  \code{\link{getIntegrationSites}},
  \code{\link{pslToRangedObject}}
}

