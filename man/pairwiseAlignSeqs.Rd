\name{pairwiseAlignSeqs}
\alias{pairwiseAlignSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pairwiseAlignSeqs(subjectSeqs = NULL, patternSeq = NULL, side = "left", qualityThreshold = 1, showStats = FALSE, bufferBases = 5, doRC = TRUE, returnUnmatched = FALSE, returnLowScored = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subjectSeqs}{
%%     ~~Describe \code{subjectSeqs} here~~
}
  \item{patternSeq}{
%%     ~~Describe \code{patternSeq} here~~
}
  \item{side}{
%%     ~~Describe \code{side} here~~
}
  \item{qualityThreshold}{
%%     ~~Describe \code{qualityThreshold} here~~
}
  \item{showStats}{
%%     ~~Describe \code{showStats} here~~
}
  \item{bufferBases}{
%%     ~~Describe \code{bufferBases} here~~
}
  \item{doRC}{
%%     ~~Describe \code{doRC} here~~
}
  \item{returnUnmatched}{
%%     ~~Describe \code{returnUnmatched} here~~
}
  \item{returnLowScored}{
%%     ~~Describe \code{returnLowScored} here~~
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
function (subjectSeqs = NULL, patternSeq = NULL, side = "left", 
    qualityThreshold = 1, showStats = FALSE, bufferBases = 5, 
    doRC = TRUE, returnUnmatched = FALSE, returnLowScored = FALSE, 
    ...) 
{
    if (is.null(subjectSeqs) | is.null(patternSeq)) {
        stop("subjectSeqs/patternSeq is empty. Please supply reads to be aligned")
    }
    if (length(patternSeq) > 1) {
        stop("More than 1 patternSeq is defined. Please only supply one pattern.")
    }
    if (is.null(names(subjectSeqs))) {
        names(subjectSeqs) <- paste("read", 1:length(subjectSeqs))
    }
    if (tolower(side) == "left") {
        badSeqs <- DNAStringSet()
        culprits <- width(subjectSeqs) < (nchar(patternSeq) + 
            bufferBases)
        if (any(culprits)) {
            badSeqs <- subjectSeqs[culprits]
            message(length(badSeqs), " sequences were removed from aligning since they were shorter than pattern getting aligned: ", 
                (nchar(patternSeq) + bufferBases), "bp")
            subjectSeqs <- subjectSeqs[!culprits]
        }
        subjectSeqs2 <- subseq(subjectSeqs, start = 1, end = (nchar(patternSeq) + 
            bufferBases))
        overFromLeft <- rep(0, length(subjectSeqs))
    }
    else if (tolower(side) == "right") {
        overFromLeft <- width(subjectSeqs) - (nchar(patternSeq) + 
            bufferBases)
        overFromLeft[overFromLeft < 1] <- 1
        subjectSeqs2 <- subseq(subjectSeqs, start = overFromLeft)
    }
    else {
        subjectSeqs2 <- subjectSeqs
        overFromLeft <- rep(0, length(subjectSeqs))
    }
    if (doRC) {
        patternSeq <- doRCtest(subjectSeqs2, patternSeq, qualityThreshold)
    }
    hits <- pairwiseAlignment(subjectSeqs2, patternSeq, type = "overlap", 
        gapOpening = -1, gapExtension = -1, ...)
    stopifnot(length(hits) == length(subjectSeqs2))
    scores <- round(score(hits))
    highscored <- scores >= round(nchar(patternSeq) * qualityThreshold) * 
        2
    unmatched <- nchar(hits) <= round(nchar(patternSeq) * 0.1)
    if (showStats & qualityThreshold != 0) {
        message("Total of ", as.numeric(table(highscored)["FALSE"]), 
            " did not have the defined pattern sequence (", patternSeq, 
            ") that passed qualityThreshold on the ", side, " side")
    }
    starts <- start(pattern(hits))
    ends <- end(pattern(hits))
    namesq <- names(subjectSeqs)
    hits <- IRanges(start = starts + overFromLeft - ifelse(side == 
        "right", 2, 0), end = ends + overFromLeft - ifelse(side == 
        "right", 2, 0), names = namesq)
    rm("scores", "subjectSeqs2", "subjectSeqs", "starts", "ends", 
        "namesq")
    if (!returnLowScored & !returnUnmatched) {
        hits <- hits[highscored]
    }
    else {
        hitstoreturn <- IRangesList(hits = hits[highscored])
        if (returnLowScored & length(hits[!highscored]) > 0) {
            hitstoreturn <- append(hitstoreturn, IRangesList(Rejected = hits[!highscored]))
        }
        if (returnUnmatched & length(hits[unmatched]) > 0) {
            hitstoreturn <- append(hitstoreturn, IRangesList(Absent = hits[unmatched]))
        }
        hits <- hitstoreturn
        rm(hitstoreturn)
    }
    cleanit <- gc()
    return(hits)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
