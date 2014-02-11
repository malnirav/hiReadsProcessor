\name{pairUpAlignments}
\alias{pairUpAlignments}
\title{Pair upalignments in a GRanges/GAlignments object}
\usage{
pairUpAlignments(psl.rd = NULL)
}
\arguments{
  \item{psl.rd}{a GRanges/GAlignments object with qNames
  ending in "}

  \item{startWithin}{upper bound limit of where the
  alignment should start within the query. Default is 3.}

  \item{alignRatioThreshold}{cuttoff for (alignment
  span/read length). Default is 0.7.}

  \item{genomicPercentIdentity}{cuttoff for
  (1-(misMatches/matches)). Default is 0.98}
}
\value{
a GRanges/GAlignments object with paired up reads
}
\description{
Given a GRanges/GAlignment object from
\code{\link{read.BAMasPSL}}, the function uses specified
filtering parameters to pair up reads where the definition
line or qName ends with "
}
\examples{
#pairUpAlignments()
#pairUpAlignments()
}
\seealso{
\code{\link{pairwiseAlignSeqs}}, \code{\link{blatSeqs}},
\code{\link{read.blast8}}, \code{\link{read.psl}},
\code{\link{getIntegrationSites}},
\code{\link{subreadAlignSeqs}}, \code{\link{GAlignments}},
\code{\link{read.BAMasPSL}}
}

