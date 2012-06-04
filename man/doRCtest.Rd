\name{doRCtest}
\alias{doRCtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doRCtest(subjectSeqs = NULL, patternSeq = NULL, qualityThreshold = 0.5, core.use = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subjectSeqs}{
%%     ~~Describe \code{subjectSeqs} here~~
}
  \item{patternSeq}{
%%     ~~Describe \code{patternSeq} here~~
}
  \item{qualityThreshold}{
%%     ~~Describe \code{qualityThreshold} here~~
}
  \item{core.use}{
%%     ~~Describe \code{core.use} here~~
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
function (subjectSeqs = NULL, patternSeq = NULL, qualityThreshold = 0.5, 
    core.use = 2) 
{
    if (is.null(subjectSeqs) | is.null(patternSeq)) {
        stop("subjectSeqs/patternSeq is empty. Please supply reads to be aligned")
    }
    if (core.use == 1) {
        registerDoSEQ()
    }
    patternSeq.rc <- as.character(reverseComplement(DNAString(patternSeq)))
    hits <- foreach(x = iter(c(patternSeq, patternSeq.rc)), .inorder = TRUE, 
        .errorhandling = "pass", .export = c("subjectSeqs", "qualityThreshold"), 
        .packages = "Biostrings") \%dopar\% {
        counts <- pmin(vcountPattern(x, subjectSeqs, max.mismatch = round(nchar(x) * 
            (1 - qualityThreshold))), 1)
        sum(counts)
    }
    if (all(hits == 0)) {
        stop("No hits found")
    }
    if (hits[[1]] < hits[[2]]) {
        message("There were less/no good hits found for pattern ", 
            patternSeq, "\nthan its reverse complement...", patternSeq.rc, 
            "\nUsing the latter to perform the searching.")
        patternSeq <- patternSeq.rc
    }
    return(patternSeq)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
