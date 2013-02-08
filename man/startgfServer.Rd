\name{startgfServer}
\alias{startgfServer}
\title{Start a gfServer instance}
\usage{
  startgfServer(seqDir = NULL, host = "localhost",
    port = 5560,
    gfServerOpts = c(repMatch = 112312, stepSize = 5, tileSize = 10))
}
\arguments{
  \item{seqDir}{absolute or relative path to the genome
  index (nib/2bit files).}

  \item{host}{name of the machine to run gfServer on.
  Default: localhost}

  \item{port}{a port number to host the gfServer with.
  Default is 5560.}

  \item{gfServerOpts}{a character vector of options to be
  passed to gfServer command on top of server defaults.
  Default: c(repMatch=112312, stepSize=5, tileSize=10). Set
  this to NULL to start gfServer with defaults.}
}
\value{
  system command status for executing gfServer command.
}
\description{
  Start a gfServer indexed reference genome to align batch
  of sequences using BLAT gfServer/gfClient protocol.
}
\examples{
#startgfServer(seqDir="/usr/local/blatSuite34/hg18.2bit",port=5560)
#stopgfServer(port=5560)
}
\seealso{
  \code{\link{stopgfServer}}, \code{\link{read.psl}},
  \code{\link{blatSeqs}}, \code{\link{read.blast8}}
}

