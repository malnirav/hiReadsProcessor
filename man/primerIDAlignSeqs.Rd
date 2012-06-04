\name{primerIDAlignSeqs}
\alias{primerIDAlignSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
primerIDAlignSeqs(subjectSeqs = NULL, patternSeq = NULL, qualityThreshold1 = 0.75, qualityThreshold2 = 0.5, doAnchored = FALSE, doRC = TRUE, returnUnmatched = FALSE, returnRejected = FALSE, showStats = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subjectSeqs}{
%%     ~~Describe \code{subjectSeqs} here~~
}
  \item{patternSeq}{
%%     ~~Describe \code{patternSeq} here~~
}
  \item{qualityThreshold1}{
%%     ~~Describe \code{qualityThreshold1} here~~
}
  \item{qualityThreshold2}{
%%     ~~Describe \code{qualityThreshold2} here~~
}
  \item{doAnchored}{
%%     ~~Describe \code{doAnchored} here~~
}
  \item{doRC}{
%%     ~~Describe \code{doRC} here~~
}
  \item{returnUnmatched}{
%%     ~~Describe \code{returnUnmatched} here~~
}
  \item{returnRejected}{
%%     ~~Describe \code{returnRejected} here~~
}
  \item{showStats}{
%%     ~~Describe \code{showStats} here~~
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
function (subjectSeqs = NULL, patternSeq = NULL, qualityThreshold1 = 0.75, 
    qualityThreshold2 = 0.5, doAnchored = FALSE, doRC = TRUE, 
    returnUnmatched = FALSE, returnRejected = FALSE, showStats = FALSE, 
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
    if (length(unlist(gregexpr("N", patternSeq))) < 4) {
        stop("There should be minimum of atleast 4 Ns in patternSeq to be considered as a primerID sequence.")
    }
    if (doRC) {
        patternSeq <- doRCtest(subjectSeqs, patternSeq)
    }
    primerIDpos <- unlist(gregexpr("N", patternSeq))
    pattern1 <- as.character(subseq(DNAString(patternSeq), 1, 
        primerIDpos[1] - 1))
    pattern2 <- as.character(subseq(DNAString(patternSeq), primerIDpos[length(primerIDpos)] + 
        1))
    pattern1.hits <- pairwiseAlignSeqs(subjectSeqs, pattern1, 
        "middle", qualityThreshold = qualityThreshold1, doRC = FALSE, 
        returnUnmatched = TRUE, ...)
    pattern2.hits <- pairwiseAlignSeqs(subjectSeqs, pattern2, 
        "middle", qualityThreshold = qualityThreshold2, doRC = FALSE, 
        ...)
    unmatched <- pattern1.hits[["Absent"]]
    pattern1.hits <- pattern1.hits[["hits"]]
    rejected1 <- setdiff(names(pattern1.hits), names(pattern2.hits))
    rejected2 <- setdiff(names(pattern2.hits), names(pattern1.hits))
    rejected <- c(pattern1.hits[names(pattern1.hits) \%in\% rejected1], 
        pattern2.hits[names(pattern2.hits) \%in\% rejected2])
    rm("rejected1", "rejected2")
    if (showStats) {
        message("Removed ", length(rejected), " read(s) for only matching one of pattern1 or pattern2")
    }
    good.rows <- intersect(names(pattern1.hits), names(pattern2.hits))
    pattern1.hits <- pattern1.hits[names(pattern1.hits) \%in\% 
        good.rows]
    pattern2.hits <- pattern2.hits[names(pattern2.hits) \%in\% 
        good.rows]
    stopifnot(identical(names(pattern1.hits), names(pattern2.hits)))
    stopifnot(identical(names(pattern1.hits), good.rows))
    badAss <- end(pattern1.hits) >= start(pattern2.hits)
    if (any(badAss)) {
        rejected <- c(rejected, pattern1.hits[badAss], pattern2.hits[badAss])
        pattern1.hits <- pattern1.hits[!badAss]
        pattern2.hits <- pattern2.hits[!badAss]
        good.rows <- good.rows[!badAss]
        message("Removed ", table(badAss)["TRUE"], " read(s) for not having primerID present between pattern1 and pattern2")
    }
    hits <- IRanges(start = start(pattern1.hits), end = end(pattern2.hits), 
        names = good.rows)
    primerIDs <- IRanges(start = end(pattern1.hits) + 1, end = start(pattern2.hits) - 
        1, names = good.rows)
    if (length(hits) == 0) {
        stop("No hits found that matched both sides of patternSeq with primerID in the middle.")
    }
    if (doAnchored) {
        message("Found ", length(primerIDs), " total primerIDs before anchored filter.")
        anchorBase.s <- substr(patternSeq, primerIDpos[1] - 1, 
            primerIDpos[1] - 1)
        anchorBase.e <- substr(patternSeq, primerIDpos[length(primerIDpos)] + 
            1, primerIDpos[length(primerIDpos)] + 1)
        anchorBase.s.i <- trimSeqs(subjectSeqs, resize(pattern1.hits, 
            width = 1, fix = "end"))
        anchorBase.e.i <- trimSeqs(subjectSeqs, resize(pattern2.hits, 
            width = 1, fix = "start"))
        rows <- anchorBase.s == as.character(anchorBase.s.i) & 
            anchorBase.e == as.character(anchorBase.e.i)
        unAnchored <- hits[!rows]
        primerIDs <- primerIDs[rows]
        hits <- hits[rows]
        message("Found ", length(primerIDs), " total primerIDs after anchored filter.")
        rm("rows", "anchorBase.s.i", "anchorBase.e.i")
    }
    rm("good.rows", "pattern1.hits", "pattern2.hits")
    cleanit <- gc()
    nSize <- 2
    badAss <- !width(primerIDs) \%in\% (length(primerIDpos) - nSize):(length(primerIDpos) + 
        nSize)
    if (any(badAss)) {
        rejectedprimerIDs <- hits[badAss]
        hits <- hits[!badAss]
        primerIDs <- primerIDs[!badAss]
        message("Removed ", table(badAss)["TRUE"], " read(s) for not having right primerID length")
    }
    hits <- IRangesList(hits = hits, primerIDs = primerIDs)
    if (exists("unAnchored")) {
        if (length(unAnchored) > 0) {
            hits <- append(hits, IRangesList(unAnchoredprimerIDs = unAnchored))
        }
    }
    if (returnUnmatched & length(unmatched) > 0) {
        hits <- append(hits, IRangesList(Absent = unmatched))
    }
    if (returnRejected) {
        if (length(rejected) > 0) {
            hits <- append(hits, IRangesList(Rejected = rejected))
        }
        if (exists("rejectedprimerIDs")) {
            if (length(rejectedprimerIDs) > 0) {
                hits <- append(hits, IRangesList(RejectedprimerIDs = rejectedprimerIDs))
            }
        }
    }
    return(hits)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
