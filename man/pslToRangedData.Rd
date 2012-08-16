\name{pslToRangedData}
\alias{pslToRangedData}
\title{Convert psl dataframe to RangedData}
\usage{
  pslToRangedData(x, useTargetAsRef = TRUE)
}
\arguments{
  \item{x}{dataframe reflecting psl format}

  \item{useTargetAsRef}{use target or query as space or the
  reference data. Default is TRUE.}
}
\value{
  a RangedData object reflecting psl file type.
}
\description{
  Convert psl dataframe to RangedData object using either
  the query or target as the reference data column.
}
\examples{
#pslToRangedData(psl)
 #pslToRangedData(psl,useTargetAsRef=FALSE)
}
\seealso{
  \code{\link{blatListedSet}}
}

