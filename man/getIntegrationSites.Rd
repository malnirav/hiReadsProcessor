\name{getIntegrationSites}
\alias{getIntegrationSites}
\title{Obtain integration sites from BLAT output.
Given a RangedData object from \code{\link{read.psl}}, the function uses specified filtering parameters to obtain integration sites and maintain sequence attrition. The function will remove any non-best scoring alignments from the object if not already filtered apriori.}
\usage{
  getIntegrationSites(psl.rd = NULL, startWithin = 3,
    alignRatioThreshold = 0.7,
    genomicpercentidentity = 0.98, correctByqStart = TRUE)
}
\arguments{
  \item{psl.rd}{a RangedData object reflecting psl format
  where tName is the spaces.}

  \item{startWithin}{upper bound limit of where the
  alignment should start within the query. Default is 3.}

  \item{alignRatioThreshold}{cuttoff for (alignment
  span/read length). Default is 0.7.}

  \item{genomicpercentidentity}{cuttoff for
  (1-(misMatches/matches)). Default is 0.98.}

  \item{correctByqStart}{use qStart to correct genomic
  position. This would account for sequencing/trimming
  errors.
  Position=ifelse(strand=="+",tStart-qStart,tEnd+qStart).
  Default is TRUE.}
}
\value{
  a RangedData object with integration sites which passed
  all filtering criteria. Each filtering parameter creates
  a new column to flag if a sequence/read passed that
  filter which follows the scheme: 'pass.FilterName'.
}
\description{
  Obtain integration sites from BLAT output. Given a
  RangedData object from \code{\link{read.psl}}, the
  function uses specified filtering parameters to obtain
  integration sites and maintain sequence attrition. The
  function will remove any non-best scoring alignments from
  the object if not already filtered apriori.
}
\examples{
\dontrun{
 getIntegrationSites(test.psl.rd)
}
}
\seealso{
  \code{\link{startgfServer}}, \code{\link{read.psl}},
  \code{\link{blatSeqs}}, \code{\link{blatListedSet}},
  \code{\link{findIntegrations}},
  \code{\link{pslToRangedData}},
  \code{\link{clusterSites}}, \code{\link{otuSites}}
}

