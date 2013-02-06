\name{pslToRangedObject}
\alias{pslToRangedObject}
\title{Convert psl dataframe to RangedData/GRanges}
\usage{
  pslToRangedObject(x, useTargetAsRef = TRUE,
    asGRanges = FALSE, isblast8 = FALSE)
}
\arguments{
  \item{x}{dataframe reflecting psl format}

  \item{useTargetAsRef}{use target or query as space or the
  reference data. Default is TRUE.}

  \item{asGRanges}{make a GRanges object instead of
  RangedData. Default is FALSE.}

  \item{isblast8}{the input dataframe blast8 format output
  from BLAT. Default is FALSE.}
}
\value{
  a RangedData/GRanges object reflecting psl file type.
}
\description{
  Convert psl dataframe to RangedData or GRanges object
  using either the query or target as the reference data
  column.
}
\examples{
#pslToRangedObject(psl)
#pslToRangedObject(psl,asGRanges=TRUE)
#pslToRangedObject(psl,useTargetAsRef=FALSE)
}
\seealso{
  \code{\link{blatListedSet}}
}

