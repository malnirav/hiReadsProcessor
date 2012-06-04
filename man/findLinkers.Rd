\name{findLinkers}
\alias{findLinkers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findLinkers(sampleInfo, showStats = FALSE, doRC = FALSE, parallel = TRUE, samplenames = NULL)
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
    primerFlag <- FALSE
    toProcess <- extractFeature(sampleInfo, feature = "LTRed")
    toProcessSamples <- sapply(toProcess, names, simplify = FALSE)
    sectors <- names(which(sapply(sapply(toProcess, length), 
        ">", 1)))
    if (length(sectors) == 0) {
        message("No LTRed information found in sampleInfo. Using primed information to find Linkers")
        toProcess <- extractFeature(sampleInfo, feature = "primed")
        toProcessSamples <- sapply(toProcess, names, simplify = FALSE)
        sectors <- names(which(sapply(sapply(toProcess, length), 
            ">", 1)))
        if (length(sectors) == 0) {
            stop("No primed information found in sampleInfo...did you run findPrimers()?")
        }
        primerFlag <- TRUE
    }
    rm(toProcess)
    gc()
    for (sector in sectors) {
        message("Processing sector ", sector)
        sampleLinkers <- toupper(extractFeature(sampleInfo, sector = sector, 
            feature = "linkersequence")[[sector]])
        if (length(sampleLinkers) == 0 | mean(nchar(sampleLinkers)) <= 
            10 | all(is.na(sampleLinkers))) {
            stop("Either Linker sequence is too short (<=10) or no Linkers found in sample information object.")
        }
        linkerIdentity <- extractFeature(sampleInfo, sector = sector, 
            feature = "linkeridentity")[[sector]]
        stopifnot(length(linkerIdentity) > 0)
        primerIded <- extractFeature(sampleInfo, sector = sector, 
            feature = "primeridinlinker")[[sector]]
        primerIded.threshold1 <- extractFeature(sampleInfo, sector = sector, 
            feature = "primeridinlinkeridentity1")[[sector]]
        primerIded.threshold2 <- extractFeature(sampleInfo, sector = sector, 
            feature = "primeridinlinkeridentity2")[[sector]]
        samplesToProcess <- toProcessSamples[[sector]]
        if (!is.null(samplenames)) {
            samplesToProcess <- samplesToProcess[samplesToProcess \%in\% 
                samplenames]
        }
        toProcess <- extractSeqs(sampleInfo, sector, samplesToProcess, 
            feature = ifelse(primerFlag, "primed", "LTRed"))[[sector]]
        message("\tFinding Linkers.")
        linkerTrimmed <- foreach(x = iter(samplesToProcess), 
            .inorder = TRUE, .errorhandling = "pass", .export = c("primerIded", 
                "toProcess", "sampleLinkers", "doRC", "linkerIdentity", 
                "primerIded.threshold1", "primerIded.threshold2", 
                "pairwiseAlignSeqs", "primerIDAlignSeqs"), .packages = "Biostrings") \%dopar\% 
            {
                if (primerIded[[x]]) {
                  primerIDAlignSeqs(toProcess[[x]], sampleLinkers[[x]], 
                    doAnchored = TRUE, returnUnmatched = TRUE, 
                    returnRejected = TRUE, doRC = doRC, qualityThreshold1 = primerIded.threshold1[[x]], 
                    qualityThreshold2 = primerIded.threshold2[[x]])
                }
                else {
                  pairwiseAlignSeqs(toProcess[[x]], sampleLinkers[[x]], 
                    "middle", qualityThreshold = linkerIdentity[[x]], 
                    returnUnmatched = TRUE, returnLowScored = TRUE, 
                    doRC = doRC)
                }
            }
        names(linkerTrimmed) <- samplesToProcess
        if (any(grepl("simpleError", linkerTrimmed))) {
            stop("Error encountered in Linker Trimming functions", 
                paste(names(linkerTrimmed[grepl("simpleError", 
                  linkerTrimmed)]), collapse = ", "))
        }
        culprits <- grep("No hits found", linkerTrimmed)
        if (length(culprits) > 0) {
            message("Following sample(s) had no hits for Linker alignment: ", 
                paste(samplesToProcess[culprits], collapse = ", "))
            samplesToProcess <- samplesToProcess[-c(culprits)]
            linkerTrimmed <- linkerTrimmed[-c(culprits)]
        }
        cleanit <- gc()
        toprint <- data.frame(Total = sapply(sapply(linkerTrimmed, 
            "[[", "hits"), length))
        counts <- sapply(toProcess, length)
        toprint$PercOfPrimedOrLTRed <- 100 * (toprint$Total/counts[rownames(toprint)])
        toprint$SampleName <- rownames(toprint)
        rownames(toprint) <- NULL
        if (showStats) {
            message("Sequence Stats after Linker alignment:")
            print(toprint)
        }
        if (mean(toprint$PercOfPrimedOrLTRed) <= 5) {
            stop("Something seems to be wrong with the linkers provided for each sample. On average <= 5\% of sequences found linker match for the entire run!!!")
        }
        message("Adding linker info back to the object")
        for (x in names(linkerTrimmed)) {
            for (y in names(linkerTrimmed[[x]])) {
                cat(".")
                if (primerFlag) {
                  LTRed.ends <- end(sampleInfo$sectors[[sector]]$samples[[x]]$primed[names(linkerTrimmed[[x]][[y]])])
                }
                else {
                  LTRed.ends <- end(sampleInfo$sectors[[sector]]$samples[[x]]$LTRed[names(linkerTrimmed[[x]][[y]])])
                }
                end(linkerTrimmed[[x]][[y]]) <- end(linkerTrimmed[[x]][[y]]) + 
                  LTRed.ends
                start(linkerTrimmed[[x]][[y]]) <- start(linkerTrimmed[[x]][[y]]) + 
                  LTRed.ends
                newAttrName <- paste(ifelse(y == "hits", "", 
                  y), "linkered", sep = "")
                sampleInfo$sectors[[sector]]$samples[[x]][[newAttrName]] <- linkerTrimmed[[x]][[y]]
                rm(LTRed.ends)
            }
        }
        rm("linkerTrimmed")
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
