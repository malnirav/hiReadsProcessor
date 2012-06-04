\name{getIntegrationSites}
\alias{getIntegrationSites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getIntegrationSites(psl.rd = NULL, startWithin = 3, alignRatioThreshold = 0.7, genomicpercentidentity = 0.98, correctByqStart = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psl.rd}{
%%     ~~Describe \code{psl.rd} here~~
}
  \item{startWithin}{
%%     ~~Describe \code{startWithin} here~~
}
  \item{alignRatioThreshold}{
%%     ~~Describe \code{alignRatioThreshold} here~~
}
  \item{genomicpercentidentity}{
%%     ~~Describe \code{genomicpercentidentity} here~~
}
  \item{correctByqStart}{
%%     ~~Describe \code{correctByqStart} here~~
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
function (psl.rd = NULL, startWithin = 3, alignRatioThreshold = 0.7, 
    genomicpercentidentity = 0.98, correctByqStart = TRUE) 
{
    stopifnot(class(psl.rd) == "RangedData" & !is.null(psl.rd) & 
        !is.null(startWithin) & !is.null(alignRatioThreshold) & 
        !is.null(genomicpercentidentity))
    if (correctByqStart) {
        psl.rd$Position <- ifelse(psl.rd$strand == "+", start(psl.rd) - 
            psl.rd$qStart, end(psl.rd) + psl.rd$qStart)
    }
    else {
        psl.rd$Position <- ifelse(psl.rd$strand == "+", start(psl.rd), 
            end(psl.rd))
    }
    if (!"score" \%in\% colnames(psl.rd)) {
        psl.rd$score <- with(psl.rd, matches - misMatches - qBaseInsert - 
            tBaseInsert)
        bestScore <- tapply(psl.rd$score, space(psl.rd), max)
        isBest <- with(psl.rd, score == bestScore[space(psl.rd)])
        psl.rd <- psl.rd[isBest, ]
        rm("isBest", "bestScore")
        gc()
    }
    message("Performing QC checks.")
    psl.rd$pass.startWithin <- psl.rd$qStart <= startWithin
    psl.rd$alignRatio <- with(psl.rd, score/qSize)
    psl.rd$pass.alignRatio <- psl.rd$alignRatio >= alignRatioThreshold
    psl.rd$percIdentity <- with(psl.rd, 1 - (misMatches/matches))
    psl.rd$pass.percIdentity <- psl.rd$percIdentity >= genomicpercentidentity
    cloneHits <- table(psl.rd$qName)
    psl.rd$isMultiHit <- as.logical(cloneHits[as.character(psl.rd$qName)] > 
        1)
    rm(cloneHits)
    gc()
    psl.rd$pass.allQC <- with(psl.rd, pass.percIdentity & pass.alignRatio & 
        pass.startWithin)
    return(psl.rd)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
