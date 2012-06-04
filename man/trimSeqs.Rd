\name{trimSeqs}
\alias{trimSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trimSeqs(dnaSet, coords, side = "middle", offBy = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dnaSet}{
%%     ~~Describe \code{dnaSet} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{side}{
%%     ~~Describe \code{side} here~~
}
  \item{offBy}{
%%     ~~Describe \code{offBy} here~~
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
function (dnaSet, coords, side = "middle", offBy = 0) 
{
    stopifnot(class(dnaSet) \%in\% c("DNAStringSet", "DNAString"))
    stopifnot(class(coords) == "IRanges")
    if (is.null(names(dnaSet)) | is.null(names(coords))) {
        stopifnot(length(dnaSet) == length(coords))
    }
    else {
        rows <- match(names(coords), names(dnaSet))
        if (any(is.na(rows))) {
            stop("Some of the names in coords are not present in dnaSet")
        }
        if (!is.ordered(rows)) {
            dnaSet <- dnaSet[rows]
            if (!identical(names(dnaSet), names(coords))) {
                stop("Names are not identical between dnaSet and coords parameter")
            }
        }
    }
    if (tolower(side) == "left") {
        test <- end(coords) + offBy > width(dnaSet) | end(coords) + 
            offBy < 1
        if (any(test)) {
            message("Following sequences were removed from trimming since their coordinates+offBy were out of sequence length: ", 
                paste(names(dnaSet)[test], collapse = ", "))
        }
        subseq(dnaSet[!test], start = end(coords[!test]) + offBy)
    }
    else if (tolower(side) == "right") {
        test <- start(coords) - offBy > width(dnaSet) | end(coords) + 
            offBy < 1
        if (any(test)) {
            message("Following sequences were removed from trimming since their coordinates+offBy were out of sequence length: ", 
                paste(names(dnaSet)[test], collapse = ", "))
        }
        subseq(dnaSet[!test], end = start(coords[!test]) - offBy)
    }
    else {
        test <- start(coords) + offBy > width(dnaSet) | end(coords) - 
            offBy > width(dnaSet) | start(coords) + offBy < 1
        if (any(test)) {
            message("Following sequences were removed from trimming since their coordinates+offBy were out of sequence length: ", 
                paste(names(dnaSet)[test], collapse = ", "))
        }
        subseq(dnaSet[!test], start = start(coords[!test]) + 
            offBy, end = end(coords[!test]) - offBy)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
