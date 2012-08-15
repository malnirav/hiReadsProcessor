\name{otuSites}
\alias{otuSites}
\title{Make OTUs of discrete positions grouped by reads}
\usage{
  otuSites(posID = NULL, readID = NULL, grouping = NULL,
    psl.rd = NULL)
}
\arguments{
  \item{posID}{a vector of discrete positions, i.e.
  Chr,strand,Position.}

  \item{readID}{a vector of read/clone names which is
  unique to each row, i.e. deflines.}

  \item{grouping}{additional vector of grouping by which to
  pool the rows (i.e. samplenames). Default is NULL.}

  \item{psl.rd}{a RangedData object returned from
  \code{\link{clusterSites}}. Default is NULL.}
}
\value{
  a data frame with posID, readID, grouping, and otuID. If
  psl.rd parameter is defined, then a RangedData object
  where object is first filtered by clusterTopHit column
  and the otuID column appended at the end.
}
\description{
  Given a group of discrete positions per read/clone, the
  function tries to yield a unique OTU ID for the
  collection based on overlap of discrete positions to
  other reads/clones by grouping. This is mainly useful
  when each readID has many posID which needs to be
  considered as one single group of sites.
}
\examples{
\dontrun{
 otuSites(posID=c('chr1-1000','chr1-1000','chr2-1000','chr2+1000','chr15-1000','chr16-1000','chr11-1000'), readID=paste('read',sample(letters,7),sep='-'), grouping=c('a','a','a','b','b','b','c'))
 otuSites(psl.rd=test.psl.rd)
}
}
\seealso{
  \code{\link{startgfServer}}, \code{\link{read.psl}},
  \code{\link{blatSeqs}}, \code{\link{blatListedSet}},
  \code{\link{findIntegrations}},
  \code{\link{pslToRangedData}},
  \code{\link{getIntegrationSites}},
  \code{\link{clusterSites}}
}

