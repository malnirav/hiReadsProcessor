\name{replicateReads}
\alias{replicateReads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replicateReads(dnaSet, counts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dnaSet}{
%%     ~~Describe \code{dnaSet} here~~
}
  \item{counts}{
%%     ~~Describe \code{counts} here~~
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
function (dnaSet, counts = NULL) 
{
    stopifnot(class(dnaSet) == "DNAStringSet")
    if (is.null(counts)) {
        if (is.null(names(dnaSet))) {
            stop("No names attribute found in dnaSet object")
        }
        counts <- as.numeric(sub(".+counts=(.+)", "\\1", names(dnaSet)))
        if (all(is.na(counts))) {
            stop("No counts=X marker found at the end of definition line or names attribute in dnaSet object")
        }
    }
    if (length(counts) == 1) {
        counts <- rep(counts, length(dnaSet))
    }
    ids <- unlist(sapply(counts, function(x) 1:x))
    deflines <- paste(rep(sub("(.+)counts=.+", "\\1", names(dnaSet)), 
        times = counts), ids, sep = "")
    dnaSet <- rep(dnaSet, times = counts)
    names(dnaSet) <- deflines
    return(dnaSet)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
