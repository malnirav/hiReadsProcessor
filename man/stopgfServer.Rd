\name{stopgfServer}
\alias{stopgfServer}
\title{Stop a gfServer instance
Take down the gfServer instance.}
\usage{
  stopgfServer(host = "localhost", port = NULL)
}
\arguments{
  \item{host}{name of the machine running gfServer.
  Default: localhost}

  \item{port}{the same number you started the gfServer
  with.}
}
\value{
  system command status for executing kill command.
}
\description{
  Stop a gfServer instance Take down the gfServer instance.
}
\examples{
\dontrun{
	stopgfServer(port=5560)
}
}
\seealso{
  \code{\link{startgfServer}}
}

