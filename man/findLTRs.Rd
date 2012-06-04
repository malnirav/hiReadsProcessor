\name{findLTRs}
\alias{findLTRs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findLTRs(sampleInfo, showStats = FALSE, doRC = FALSE, parallel = TRUE, samplenames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfo}{
%%     ~~Describe \code{sampleInfo} here~~
}
  \item{showStats}{
%%     ~~Describe \code{showStats} here~~
}
  \item{doRC}{
%%     ~~Describe \code{doRC} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{samplenames}{
%%     ~~Describe \code{samplenames} here~~
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
function (sampleInfo, showStats = FALSE, doRC = FALSE, parallel = TRUE, 
    samplenames = NULL) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    if (!parallel) {
        registerDoSEQ()
    }
    primed <- extractFeature(sampleInfo, feature = "primed")
    samplesprimed <- sapply(primed, names, simplify = FALSE)
    sectorsPrimed <- names(which(sapply(sapply(primed, length), 
        ">", 1)))
    rm(primed)
    gc()
    if (length(sectorsPrimed) == 0) {
        stop("No primed information found in sampleInfo...did you run findPrimers()?")
    }
    for (sector in sectorsPrimed) {
        message("Processing sector ", sector)
        sampleLTRbits <- toupper(extractFeature(sampleInfo, sector = sector, 
            feature = "ltrbitsequence")[[sector]])
        if (length(sampleLTRbits) == 0 | mean(nchar(sampleLTRbits)) <= 
            1 | all(is.na(sampleLTRbits))) {
            stop("Either LTR bit sequence is too short (<=1) or no LTR bits found in sample information object.")
        }
        samplesToProcess <- samplesprimed[[sector]]
        if (!is.null(samplenames)) {
            samplesToProcess <- samplesToProcess[samplesToProcess \%in\% 
                samplenames]
        }
        primerTrimmed <- extractSeqs(sampleInfo, sector, samplesToProcess, 
            feature = "primed")[[sector]]
        message("\tFinding LTR bits.")
        ltrBitIdentity <- extractFeature(sampleInfo, sector = sector, 
            feature = "ltrbitidentity")[[sector]]
        ltrTrimmed <- foreach(x = iter(samplesToProcess), .inorder = TRUE, 
            .errorhandling = "pass", .export = c("primerTrimmed", 
                "sampleLTRbits", "ltrBitIdentity", "doRC", "pairwiseAlignSeqs"), 
            .packages = "Biostrings") \%dopar\% {
            pairwiseAlignSeqs(primerTrimmed[[x]], sampleLTRbits[[x]], 
                "left", qualityThreshold = ltrBitIdentity[[x]], 
                doRC = doRC)
        }
        names(ltrTrimmed) <- samplesToProcess
        if (any(grepl("simpleError", ltrTrimmed))) {
            stop("Error encountered in LTR Trimming function", 
                paste(names(ltrTrimmed[grepl("simpleError", ltrTrimmed)]), 
                  collapse = ", "))
        }
        culprits <- grep("No hits found", ltrTrimmed)
        if (length(culprits) > 0) {
            message("Following sample(s) had no hits for LTR bit alignment: ", 
                paste(samplesToProcess[culprits], collapse = ", "))
            samplesToProcess <- samplesToProcess[-c(culprits)]
            ltrTrimmed <- ltrTrimmed[-c(culprits)]
        }
        cleanit <- gc()
        toprint <- as.data.frame(sapply(ltrTrimmed, length))
        names(toprint) <- "Total"
        counts <- sapply(primerTrimmed, length)
        toprint$PercOfPrimed <- 100 * (toprint$Total/counts[rownames(toprint)])
        toprint$SampleName <- rownames(toprint)
        rownames(toprint) <- NULL
        if (showStats) {
            message("Sequence Stats after LTR bit alignment:")
            print(toprint)
        }
        if (mean(toprint$PercOfPrimed) <= 5) {
            stop("Something seems to be wrong with the LTRs provided for each sample. On average <= 5\% of sequences found LTR match for the entire run!!!")
        }
        message("Adding LTR info back to the object")
        for (x in names(ltrTrimmed)) {
            cat(".")
            if (length(ltrTrimmed[[x]]) > 0) {
                primed.end <- end(sampleInfo$sectors[[sector]]$samples[[x]]$primed[names(ltrTrimmed[[x]])])
                end(ltrTrimmed[[x]]) <- end(ltrTrimmed[[x]]) + 
                  primed.end
                start(ltrTrimmed[[x]]) <- start(ltrTrimmed[[x]]) + 
                  primed.end
                sampleInfo$sectors[[sector]]$samples[[x]]$LTRed <- ltrTrimmed[[x]]
                rm(primed.end)
            }
        }
        rm("ltrTrimmed", "primerTrimmed")
        cleanit <- gc()
    }
    sampleInfo$callHistory <- append(sampleInfo$callHistory, 
        match.call())
    return(sampleInfo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
