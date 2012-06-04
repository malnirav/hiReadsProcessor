\name{pslToRangedData}
\alias{pslToRangedData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pslToRangedData(x, useTargetAsRef = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{useTargetAsRef}{
%%     ~~Describe \code{useTargetAsRef} here~~
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
function (x, useTargetAsRef = TRUE) 
{
    if (useTargetAsRef) {
        metadataCols <- c(grep("tName|tStart|tEnd|strand", names(x), 
            invert = TRUE, value = TRUE), "tStarts")
        RangedData(space = x$tName, IRanges(start = x$tStart, 
            end = x$tEnd), strand = x$strand, x[, metadataCols])
    }
    else {
        metadataCols <- c(grep("qName|qStart|qEnd|strand", names(x), 
            invert = TRUE, value = TRUE), "qStarts")
        RangedData(space = x$qName, IRanges(start = x$qStart, 
            end = x$qEnd), strand = x$strand, x[, metadataCols])
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
