\name{read.psl}
\alias{read.psl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.psl(pslFile = NULL, bestScoring = TRUE, asRangedData = FALSE, removeFile = TRUE, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pslFile}{
%%     ~~Describe \code{pslFile} here~~
}
  \item{bestScoring}{
%%     ~~Describe \code{bestScoring} here~~
}
  \item{asRangedData}{
%%     ~~Describe \code{asRangedData} here~~
}
  \item{removeFile}{
%%     ~~Describe \code{removeFile} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
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
function (pslFile = NULL, bestScoring = TRUE, asRangedData = FALSE, 
    removeFile = TRUE, parallel = FALSE) 
{
    if (is.null(pslFile)) {
        stop("pslFile parameter empty. Please supply a filename to be read.")
    }
    if (any(grepl("\\*", pslFile))) {
        pslFile <- list.files(path = dirname(pslFile), pattern = basename(pslFile), 
            full.names = TRUE)
    }
    if (length(pslFile) == 0) {
        stop("No file(s) found with given paramter in pslFile:", 
            pslFile)
    }
    if (!parallel) {
        registerDoSEQ()
    }
    cols <- c("matches", "misMatches", "repMatches", "nCount", 
        "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", 
        "strand", "qName", "qSize", "qStart", "qEnd", "tName", 
        "tSize", "tStart", "tEnd", "blockCount", "blockSizes", 
        "qStarts", "tStarts")
    cols.class <- c(rep("numeric", 8), rep("character", 2), rep("numeric", 
        3), "character", rep("numeric", 4), rep("character", 
        3))
    hits <- foreach(x = iter(pslFile), .inorder = FALSE, .export = c("cols", 
        "cols.class", "bestScoring", "removeFile")) \%dopar\% {
        message(x)
        hits.temp <- read.delim(x, header = FALSE, col.names = cols, 
            stringsAsFactors = FALSE, colClasses = cols.class)
        if (removeFile) {
            system(paste("rm", x))
        }
        if (bestScoring) {
            hits.temp$score <- with(hits.temp, matches - misMatches - 
                qBaseInsert - tBaseInsert)
            hits.temp <- hits.temp[with(hits.temp, order(qName, 
                -score)), ]
            bestScore <- with(hits.temp, tapply(score, qName, 
                max))
            isBest <- with(hits.temp, score == bestScore[qName])
            hits.temp <- hits.temp[isBest, ]
            rm("isBest", "bestScore")
        }
        hits.temp
    }
    hits <- unique(do.call(rbind, hits))
    if (nrow(hits) == 0) {
        stop("No hits found")
    }
    message("Ordering by qName and cherry picking!")
    hits <- orderBy(~qName, hits)
    if (bestScoring) {
        hits$score <- with(hits, matches - misMatches - qBaseInsert - 
            tBaseInsert)
        hits <- hits[with(hits, order(qName, -score)), ]
        bestScore <- with(hits, tapply(score, qName, max))
        isBest <- with(hits, score == bestScore[qName])
        hits <- hits[isBest, ]
        rm("isBest", "bestScore")
    }
    if (asRangedData) {
        hits <- pslToRangedData(hits, useTargetAsRef = TRUE)
    }
    return(hits)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
