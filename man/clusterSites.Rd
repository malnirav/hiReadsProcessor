\name{clusterSites}
\alias{clusterSites}
\title{Cluster values within a window based on their frequency given discrete factors}
\usage{
  clusterSites(posID = NULL, value = NULL, grouping = NULL,
    psl.rd = NULL, weight = NULL, windowSize = 5,
    byQuartile = FALSE, quartile = 0.7, parallel = TRUE)
}
\arguments{
  \item{posID}{a vector of groupings for the value
  parameter (i.e. Chr,strand). Required if psl.rd parameter
  is not defined.}

  \item{value}{a vector of integer with values that needs
  to corrected/clustered (i.e. Positions). Required if
  psl.rd parameter is not defined.}

  \item{grouping}{additional vector of grouping by which to
  pool the rows (i.e. samplenames). Default is NULL.}

  \item{psl.rd}{a RangedData object returned from
  \code{\link{getIntegrationSites}}. Default is NULL.}

  \item{weight}{a numeric vector of weights to use when
  calculating frequency of value by posID and grouping if
  specified. Default is NULL.}

  \item{windowSize}{size of window for which values should
  be corrected/clustered. Default is 5.}

  \item{byQuartile}{flag denoting whether quartile based
  technique should be employed. See notes for details.
  Default is TRUE.}

  \item{quartile}{if byQuartile=TRUE, then the quartile
  which serves as the threshold. Default is 0.70.}

  \item{parallel}{use parallel backend to perform
  calculation with \code{\link{foreach}}. Defaults to TRUE.
  If no parallel backend is registered, then a serial
  version of foreach is ran using
  \code{\link{registerDoSEQ()}}. Process is split by the
  grouping the column.}
}
\value{
  a data frame with clusteredValues and frequency shown
  alongside with the original input. If psl.rd parameter is
  defined then a RangedData object is returned with three
  new columns appended at the end: clusteredPosition,
  clonecount, and clusterTopHit (a representative for a
  given cluster chosen by best scoring hit!).
}
\description{
  Given a group of discrete factors (i.e. position ids) and
  integer values, the function tries of correct/cluster the
  integer values based on their frequency in a defined
  windowsize.
}
\note{
  The algorithm for clustering when byQuartile=TRUE is as
  follows: for all values in each grouping, get a
  distribution and test if their frequency is >= quartile
  threshold. For values below the quartile threshold, test
  if any values overlap with the ones that passed the
  threshold and is within the defined windowSize. If there
  is a match, then merge with higher value, else leave it
  as is. This is only useful if the distribution is wide
  and polynodal. When byQuartile=FALSE, for each group the
  values within the defined window are merged with the next
  highest frequently occuring value, if freuquencies are
  tied then lowest value is used to represent the cluster.
}
\examples{
#clusterSites(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), value=c(rep(1000,2),1004,1000,1000,1000,1000), grouping=c('a','a','a','b','b','b','c'))
 #clusterSites(grouping=test.psl.rd$grouping, psl.rd=test.psl.rd)
}
\seealso{
  \code{\link{startgfServer}}, \code{\link{read.psl}},
  \code{\link{blatSeqs}}, \code{\link{blatListedSet}},
  \code{\link{findIntegrations}},
  \code{\link{pslToRangedData}},
  \code{\link{getIntegrationSites}},
  \code{\link{otuSites}}, \code{\link{read.blast8}}
}

