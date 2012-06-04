\name{blatListedSet}
\alias{blatListedSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
blatListedSet(dnaSetList = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dnaSetList}{
%%     ~~Describe \code{dnaSetList} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (dnaSetList = NULL, ...) 
{
    if (is.null(dnaSetList)) {
        stop("dnaSetList is empty. Please supply a listed DNAStringSet returned from extractSeqs() to be aligned against a reference")
    }
    sapply(names(dnaSetList), function(x) {
        do.call(RangedDataList, sapply(names(dnaSetList[[x]]), 
            function(y) {
                if (length(dnaSetList[[x]][[y]]) > 0) {
                  message("BLATing ", y)
                  filenames <- blatSeqs(query = dnaSetList[[x]][[y]], 
                    ...)
                  read.psl(filenames)
                }
            }))
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
