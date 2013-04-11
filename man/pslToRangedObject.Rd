\name{pslToRangedObject}
\alias{pslToRangedObject}
\title{Convert psl dataframe to RangedData/GRanges}
\usage{
  pslToRangedObject(x, useTargetAsRef = TRUE,
    asGRanges = TRUE, isblast8 = FALSE)
}
\arguments{
  \item{x}{dataframe reflecting psl format}

  \item{useTargetAsRef}{use target(tName) or query(qName)
  as the chromosome or the reference data. Default is
  TRUE.}

  \item{asGRanges}{make a GRanges object instead of
  RangedData. Default is TRUE.}

  \item{isblast8}{the input dataframe blast8 format output
  from BLAT. Default is FALSE.}
}
\value{
  a GRanges object reflecting psl file type. If
  asGRanges=FALSE, then RangedData object is returned.
}
\description{
  Convert psl dataframe to RangedData or GRanges object
  using either the query or target as the reference data
  column.
}
\examples{
#pslToRangedObject(psl)
#pslToRangedObject(psl, asGRanges=FALSE)
#pslToRangedObject(psl, useTargetAsRef=FALSE)
}
\seealso{
  \code{\link{read.psl}}, \code{\link{read.blast8}},
  \code{\link{blatListedSet}}
}

