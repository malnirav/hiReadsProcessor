\name{vpairwiseAlignSeqs}
\alias{vpairwiseAlignSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vpairwiseAlignSeqs(subjectSeqs = NULL, patternSeq = NULL, side = "left", qualityThreshold = 1, showStats = FALSE, bufferBases = 5, doRC = TRUE, ...)
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
    doRC = TRUE, ...) 
{
    if (is.null(subjectSeqs) | is.null(patternSeq)) {
        stop("subjectSeqs/patternSeq is empty. Please supply reads to be aligned")
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
    hits <- vmatchPattern(patternSeq, subjectSeqs2, max.mismatch = round(nchar(patternSeq) * 
        (1 - qualityThreshold)), ...)
    hits <- unlist(hits, recursive = TRUE, use.names = TRUE)
    counts <- Rle(names(hits))
    if (any(runLength(counts) > 1)) {
        reduced <- reduce(GRanges(seqnames = names(hits), IRanges(start = start(hits), 
            end = end(hits))))
        counts <- seqnames(reduced)
        hits <- ranges(reduced)
        names(hits) <- as.character(seqnames(reduced))
        rm(reduced)
        if (any(runLength(counts) > 1)) {
            message(paste("More than 1 pattern (", patternSeq, 
                ") match found for:", paste(runValue(counts)[runLength(counts) > 
                  1], collapse = ",")))
            message("\nUsing the latter occuring hit as the dominant for each read.")
            toremove <- c()
            for (culprits in as.character(runValue(counts)[runLength(counts) > 
                1])) {
                rows <- which(names(hits) \%in\% culprits)
                toremove <- c(toremove, rows[1:length(rows) - 
                  1])
            }
            hits <- hits[-toremove]
            counts <- Rle(names(hits))
            if (any(runLength(counts) > 1)) {
                stop(paste("More than 1 pattern unresolved (", 
                  patternSeq, ") match found for:", paste(runValue(counts)[runLength(counts) > 
                    1], collapse = ",")))
            }
        }
    }
    good.row <- names(subjectSeqs2) \%in\% names(hits)
    if (showStats) {
        message("Total of ", as.numeric(table(good.row)["FALSE"]), 
            " did not have the defined pattern sequence (", patternSeq, 
            ") that passed qualityThreshold on the ", side, " side")
    }
    starts <- start(hits)
    ends <- end(hits)
    namesq <- names(hits)
    rm("hits", "subjectSeqs2")
    hits <- IRanges(start = starts + overFromLeft[good.row] - 
        ifelse(side == "right", 2, 0), end = ends + overFromLeft[good.row] - 
        ifelse(side == "right", 2, 0), names = namesq)
    cleanit <- gc()
    return(hits)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
