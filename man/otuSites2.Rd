\name{otuSites2}
\alias{otuSites2}
\title{Make OTUs of genomic positions grouped by reads}
\usage{
  otuSites2(position = NULL, chromosome = NULL,
    readID = NULL, strand = "*", grouping = NULL,
    psl.rd = NULL, parallel = TRUE)
}
\arguments{
  \item{position}{a vector of integer locations/positions,
  i.e. genomic location}

  \item{chromosome}{a character vector of seqnames or
  chromosomes, i.e. tname/chromosomes}

  \item{readID}{a vector of read/clone names which is
  unique to each row, i.e. deflines.}

  \item{strand}{a character vector of strand: +,-, or *
  (default)}

  \item{grouping}{additional vector of grouping by which to
  pool the rows (i.e. samplenames). Default is NULL.}

  \item{psl.rd}{a RangedData object returned from
  \code{\link{clusterSites}}. Default is NULL.}
}
\value{
  a data frame with position, chromosome, strand, readID,
  grouping, and otuID. If psl.rd parameter is defined, then
  a RangedData object where object is first filtered by
  clusterTopHit column and the otuID column appended at the
  end.
}
\description{
  Given a group of genomic positions per read/clone, the
  function tries to yield a unique OTU ID for the
  collection based on overlap of locations to other
  reads/clones by grouping. This is mainly useful when each
  readID has many locations which needs to be considered as
  one single group of sites.
}
\examples{
#otuSites(position=sample(1:100,7),chromosome=paste0('chr',sample(1:3,7,replace=T)), readID=paste('read',sample(letters,7),sep='-'), grouping=c('a','a','a','b','b','b','c'))
 #otuSites(psl.rd=test.psl.rd)
}
\seealso{
  \code{\link{startgfServer}}, \code{\link{read.psl}},
  \code{\link{blatSeqs}}, \code{\link{blatListedSet}},
  \code{\link{findIntegrations}},
  \code{\link{pslToRangedObject}},
  \code{\link{getIntegrationSites}},
  \code{\link{clusterSites}}, \code{\link{read.blast8}}
}

