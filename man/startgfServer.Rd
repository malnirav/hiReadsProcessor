\name{startgfServer}
\alias{startgfServer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
startgfServer(host = "localhost", port = 5560, seqDir = NULL, gfServerOpts = c(repMatch = 112312, tileSize = 10), waitTime = 120)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{host}{
%%     ~~Describe \code{host} here~~
}
  \item{port}{
%%     ~~Describe \code{port} here~~
}
  \item{seqDir}{
%%     ~~Describe \code{seqDir} here~~
}
  \item{gfServerOpts}{
%%     ~~Describe \code{gfServerOpts} here~~
}
  \item{waitTime}{
%%     ~~Describe \code{waitTime} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (host = "localhost", port = 5560, seqDir = NULL, gfServerOpts = c(repMatch = 112312, 
    tileSize = 10), waitTime = 120) 
{
    if (is.null(seqDir)) {
        stop("Please define the path of nib/2bit files containing the indexed reference sequence(s)")
    }
    cmd <- paste("gfServer start", host, port, paste(paste("-", 
        names(gfServerOpts), sep = ""), gfServerOpts, collapse = " ", 
        sep = "="), seqDir, "&")
    message(cmd)
    system(cmd)
    message("Loading BLAT server...")
    system(paste("sleep", waitTime))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
