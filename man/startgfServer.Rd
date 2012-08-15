\name{startgfServer}
\alias{startgfServer}
\title{Start a gfServer instance
Start a gfServer indexed reference genome to align batch of sequences using BLAT gfServer/gfClient protocol.}
\usage{
  startgfServer(host = "localhost", port = 5560,
    seqDir = NULL,
    gfServerOpts = c(repMatch = 112312, stepSize = 5, tileSize = 10),
    waitTime = 150)
}
\arguments{
  \item{host}{name of the machine to run gfServer on.
  Default: localhost}

  \item{port}{a port number to host the gfServer with.
  Default is 5560.}

  \item{seqDir}{the path of nib/2bit files relative to the
  current working dir. Default is root.}

  \item{gfServerOpts}{a character vector of options to be
  passed to gfServer command on top of server defaults.
  Default: c(repMatch=112312, stepSize=5, tileSize=10)}

  \item{waitTime}{number of seconds to wait for gfServer to
  load. Default is 150. Increase this number for larger
  genomes, but for humans it should be suffice!}
}
\value{
  system command status for executing gfServer command.
}
\description{
  Start a gfServer instance Start a gfServer indexed
  reference genome to align batch of sequences using BLAT
  gfServer/gfClient protocol.
}
\examples{
\dontrun{
 startgfServer(port=5560,seqDir="/usr/local/blatSuite34/hg18.2bit")
 stopgfServer(port=5560)
}
}
\seealso{
  \code{\link{stopgfServer}}, \code{\link{read.psl}},
  \code{\link{blatSeqs}}
}

