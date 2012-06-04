\name{findPrimers}
\alias{findPrimers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findPrimers(sampleInfo, alignWay = "slow", showStats = FALSE, doRC = FALSE, parallel = TRUE, samplenames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfo}{
%%     ~~Describe \code{sampleInfo} here~~
}
  \item{alignWay}{
%%     ~~Describe \code{alignWay} here~~
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
function (sampleInfo, alignWay = "slow", showStats = FALSE, doRC = FALSE, 
    parallel = TRUE, samplenames = NULL) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    if (!parallel) {
        registerDoSEQ()
    }
    decoded <- extractFeature(sampleInfo, feature = "decoded")
    samplesDecoded <- sapply(decoded, names, simplify = FALSE)
    sectorsDecoded <- names(which(sapply(sapply(decoded, length), 
        ">", 1)))
    rm(decoded)
    gc()
    if (length(sectorsDecoded) == 0) {
        stop("No decoded information found in sampleInfo...did you run decodeByBarcode()?")
    }
    for (sector in sectorsDecoded) {
        message("Processing sector ", sector)
        ltrPrimers <- toupper(extractFeature(sampleInfo, sector = sector, 
            feature = "primerltrsequence")[[sector]])
        if (length(ltrPrimers) == 0 | mean(nchar(ltrPrimers)) <= 
            5 | all(is.na(ltrPrimers))) {
            stop("Either the primer size is too short (<=5) or no primers are found in sample information object.")
        }
        samplesToProcess <- samplesDecoded[[sector]]
        if (!is.null(samplenames)) {
            samplesToProcess <- samplesToProcess[samplesToProcess \%in\% 
                samplenames]
        }
        decoded <- extractSeqs(sampleInfo, sector, samplesToProcess, 
            feature = "decoded")[[sector]]
        message("\tFinding Primers.")
        primerIdentity <- extractFeature(sampleInfo, sector = sector, 
            feature = "primerltridentity")[[sector]]
        stopifnot(length(primerIdentity) > 0)
        primerTrimmed <- foreach(x = iter(samplesToProcess), 
            .inorder = TRUE, .errorhandling = "pass", .export = c("ltrPrimers", 
                "primerIdentity", "doRC", "alignWay", "vpairwiseAlignSeqs", 
                "pairwiseAlignSeqs", "decoded"), .packages = "Biostrings") \%dopar\% 
            {
                switch(alignWay, fast = vpairwiseAlignSeqs(decoded[[x]], 
                  ltrPrimers[[x]], "left", qualityThreshold = (primerIdentity[[x]] - 
                    0.05), doRC = doRC), slow = pairwiseAlignSeqs(decoded[[x]], 
                  ltrPrimers[[x]], "left", qualityThreshold = (primerIdentity[[x]]), 
                  doRC = doRC))
            }
        names(primerTrimmed) <- samplesToProcess
        if (any(grepl("simpleError", primerTrimmed))) {
            stop("Error encountered in LTR Trimming function", 
                paste(names(primerTrimmed[grepl("simpleError", 
                  primerTrimmed)]), collapse = ", "))
        }
        culprits <- grep("No hits found", primerTrimmed)
        if (length(culprits) > 0) {
            message("Following sample(s) had no hits for primer alignment: ", 
                paste(samplesToProcess[culprits], collapse = ", "))
            samplesToProcess <- samplesToProcess[-c(culprits)]
            primerTrimmed <- primerTrimmed[-c(culprits)]
        }
        cleanit <- gc()
        toprint <- as.data.frame(sapply(primerTrimmed, length))
        names(toprint) <- "Total"
        counts <- sapply(decoded, length)
        toprint$PercOfDecoded <- 100 * (toprint$Total/counts[rownames(toprint)])
        toprint$SampleName <- rownames(toprint)
        rownames(toprint) <- NULL
        if (showStats) {
            message("Sequence Stats after primer alignment:")
            print(toprint)
        }
        if (mean(toprint$PercOfDecoded) <= 5) {
            stop("Something seems to be wrong with the primers provided for each sample. On average <= 5\% of sequences found primer match for the entire run!!!")
        }
        message("Adding primer info back to the object")
        sampleInfo <- addFeature(sampleInfo, sector, names(primerTrimmed), 
            feature = "primed", value = primerTrimmed)
        rm(primerTrimmed)
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
