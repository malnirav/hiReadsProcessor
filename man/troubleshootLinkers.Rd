\name{troubleshootLinkers}
\alias{troubleshootLinkers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
troubleshootLinkers(sampleInfo, qualityThreshold = 0.55, qualityThreshold1 = 0.75, qualityThreshold2 = 0.5, doRC = TRUE, parallel = TRUE, samplenames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfo}{
%%     ~~Describe \code{sampleInfo} here~~
}
  \item{qualityThreshold}{
%%     ~~Describe \code{qualityThreshold} here~~
}
  \item{qualityThreshold1}{
%%     ~~Describe \code{qualityThreshold1} here~~
}
  \item{qualityThreshold2}{
%%     ~~Describe \code{qualityThreshold2} here~~
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
function (sampleInfo, qualityThreshold = 0.55, qualityThreshold1 = 0.75, 
    qualityThreshold2 = 0.5, doRC = TRUE, parallel = TRUE, samplenames = NULL) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    if (!parallel) {
        registerDoSEQ()
    }
    primerFlag <- FALSE
    toProcess <- extractFeature(sampleInfo, feature = "LTRed")
    toProcessSamples <- sapply(toProcess, names, simplify = FALSE)
    sectors <- names(toProcess)
    if (length(sectors) == 0) {
        message("No LTRed information found in sampleInfo. Using primed information to find Linkers")
        toProcess <- extractFeature(sampleInfo, feature = "primed")
        toProcessSamples <- sapply(toProcess, names, simplify = FALSE)
        sectors <- names(toProcess)
        if (length(sectors) == 0) {
            stop("No primed information found in sampleInfo...did you run findPrimers()?")
        }
        primerFlag <- TRUE
    }
    rm(toProcess)
    gc()
    results <- data.frame()
    for (sector in sectors) {
        message("Processing sector ", sector)
        sampleLinkers <- toupper(extractFeature(sampleInfo, sector = sector, 
            feature = "linkersequence")[[sector]])
        if (length(sampleLinkers) == 0 | mean(nchar(sampleLinkers)) <= 
            10 | all(is.na(sampleLinkers))) {
            stop("Either Linker sequence is too short (<=10) or no Linkers found in sample information object.")
        }
        samplesToProcess <- toProcessSamples[[sector]]
        if (!is.null(samplenames)) {
            samplesToProcess <- samplesToProcess[samplesToProcess \%in\% 
                samplenames]
        }
        toProcess.seqs <- extractSeqs(sampleInfo, sector, samplesToProcess, 
            feature = ifelse(primerFlag, "primed", "LTRed"))[[sector]]
        totalSeqs <- sapply(toProcess.seqs, length)
        message("\tFinding Linkers.")
        for (linkerSeq in unique(as.character(sampleLinkers))) {
            message("Checking ", linkerSeq)
            linkerTrimmed <- foreach(x = iter(samplesToProcess), 
                .inorder = TRUE, .errorhandling = "pass", .export = c("linkerSeq", 
                  "toProcess.seqs", "doRC", "qualityThreshold", 
                  "qualityThreshold1", "qualityThreshold2", "pairwiseAlignSeqs", 
                  "primerIDAlignSeqs"), .packages = "Biostrings") \%dopar\% 
                {
                  if (length(unlist(gregexpr("N", linkerSeq))) > 
                    3) {
                    length(primerIDAlignSeqs(toProcess.seqs[[x]], 
                      linkerSeq, doRC = doRC, qualityThreshold1 = qualityThreshold1, 
                      qualityThreshold2 = qualityThreshold2)$hits)
                  }
                  else {
                    length(pairwiseAlignSeqs(toProcess.seqs[[x]], 
                      linkerSeq, "middle", qualityThreshold = qualityThreshold, 
                      doRC = doRC))
                  }
                }
            names(linkerTrimmed) <- samplesToProcess
            culprits <- grep("No hits found", linkerTrimmed)
            if (length(culprits) > 0) {
                linkerTrimmed[culprits] <- 0
            }
            results <- rbind(results, data.frame(linkerSeq = linkerSeq, 
                samplename = names(linkerTrimmed), linkerhits = as.numeric(unlist(linkerTrimmed)), 
                PercentOfTotal = as.numeric(unlist(linkerTrimmed))/totalSeqs[names(linkerTrimmed)], 
                stringsAsFactors = FALSE))
            cleanit <- gc()
        }
    }
    sampleLinkers <- extractFeature(sampleInfo, feature = "linkersequence")
    names(sampleLinkers) <- NULL
    sampleLinkers <- unlist(sampleLinkers)
    linkersample <- as.data.frame(sampleLinkers)
    linkersample <- tapply(rownames(linkersample), linkersample$sampleLinkers, 
        paste, collapse = ",")
    results$CorrectLinker <- with(results, sampleLinkers[as.character(samplename)] == 
        as.character(linkerSeq))
    results$CorrectSample <- with(results, linkersample[as.character(linkerSeq)])
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
