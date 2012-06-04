\name{findIntegrations}
\alias{findIntegrations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findIntegrations(sampleInfo, seqType = NULL, port = 5560, host = "localhost", genomeIndices = c(hg18 = "/usr/local/blatSuite34/hg18.2bit", mm8 = "/usr/local/blatSuite34/mm8.2bit"), parallel = TRUE, samplenames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfo}{
%%     ~~Describe \code{sampleInfo} here~~
}
  \item{seqType}{
%%     ~~Describe \code{seqType} here~~
}
  \item{port}{
%%     ~~Describe \code{port} here~~
}
  \item{host}{
%%     ~~Describe \code{host} here~~
}
  \item{genomeIndices}{
%%     ~~Describe \code{genomeIndices} here~~
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
function (sampleInfo, seqType = NULL, port = 5560, host = "localhost", 
    genomeIndices = c(hg18 = "/usr/local/blatSuite34/hg18.2bit", 
        mm8 = "/usr/local/blatSuite34/mm8.2bit"), parallel = TRUE, 
    samplenames = NULL) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    if (!parallel) {
        registerDoSEQ()
    }
    message("Checking for linkered reads.")
    linkered <- extractFeature(sampleInfo, feature = "linkered")
    sampleslinkered <- sapply(linkered, names, simplify = FALSE)
    sectorslinkered <- names(which(sapply(sapply(linkered, length), 
        ">", 1)))
    rm(linkered)
    cleanit <- gc()
    if (length(sectorslinkered) == 0) {
        stop("No linkered information found in sampleInfo...did you run findLinkers()?")
    }
    message("Creating settings hashes for blatting and processing.")
    samplesToProcess <- unlist(sampleslinkered, use.names = F)
    if (!is.null(samplenames)) {
        samplesToProcess <- samplesToProcess[samplesToProcess \%in\% 
            samplenames]
    }
    for (setting in c("restrictionenzyme", "freeze", "startwithin", 
        "alignratiothreshold", "genomicpercentidentity", "clustersiteswithin", 
        "keepmultihits")) {
        setter <- extractFeature(sampleInfo, samplename = samplesToProcess, 
            feature = setting)
        names(setter) <- NULL
        setter <- unlist(setter)
        assign(setting, setter)
    }
    pslFiles <- c()
    for (f in unique(freeze)) {
        message("Blatting to: ", f)
        if (is.na(genomeIndices[f])) {
            stop("No genome index file defined for freeze: ", 
                f)
        }
        startgfServer(host = host, port = port, seqDir = genomeIndices[f], 
            waitTime = 120)
        message("Getting sequences to BLAT")
        if (is.null(seqType)) {
            wanted <- names(restrictionenzyme[!grepl("FRAG", 
                restrictionenzyme, ignore.case = TRUE)])
            wanted <- wanted[wanted \%in\% names(freeze[freeze == 
                f])]
            seqs <- extractSeqs(sampleInfo, samplename = wanted, 
                feature = "genomic", minReadLength = 5)
            if (any(as.numeric(sapply(seqs, length)) > 0)) {
                write.listedDNAStringSet(seqs, filePrefix = paste("processed", 
                  f, sep = ""))
            }
            wanted <- names(restrictionenzyme[grepl("FRAG", restrictionenzyme, 
                ignore.case = TRUE)])
            wanted <- wanted[wanted \%in\% names(freeze[freeze == 
                f])]
            seqs <- extractSeqs(sampleInfo, samplename = wanted, 
                feature = "genomicLinkered", minReadLength = 5)
            if (any(as.numeric(sapply(seqs, length)) > 0)) {
                write.listedDNAStringSet(seqs, filePrefix = paste("processed", 
                  f, sep = ""))
            }
        }
        else {
            seqs <- extractSeqs(sampleInfo, samplename = names(freeze[freeze == 
                f]), feature = seqType, minReadLength = 5)
            if (any(as.numeric(sapply(seqs, length)) > 0)) {
                write.listedDNAStringSet(seqs, filePrefix = paste("processed", 
                  f, sep = ""))
            }
        }
        pslFile <- blatSeqs(query = paste("processed", f, ".*.fa$", 
            sep = ""), subject = genomeIndices[[f]], standaloneBlat = FALSE, 
            host = host, port = port, parallel = parallel, gzipResults = TRUE)
        message("Kill gfServer.")
        stopgfServer(port = port)
        message("Cleaning!")
        pslFiles <- c(pslFiles, pslFile)
        cleanit <- gc()
        system(paste("rm", paste("processed", f, ".*.fa", sep = "")))
    }
    message("Reading PSL files.")
    psl <- read.psl(pslFiles, bestScoring = TRUE, asRangedData = TRUE, 
        removeFile = TRUE, parallel = FALSE)
    cleanit <- gc()
    psl$setname <- sub("^(.+)-(.+)$", "\\1", psl$qName)
    psl <- split(psl, psl$setname)
    psl.hits <- foreach(x = iter(names(psl)), .inorder = TRUE, 
        .export = c("psl", "sampleInfo", "startwithin", "alignratiothreshold", 
            "genomicpercentidentity", "clustersiteswithin", "keepmultihits")) \%dopar\% 
        {
            message("Processing ", x)
            psl.x <- getIntegrationSites(psl[[x]], startWithin = startwithin[[x]], 
                alignRatioThreshold = alignratiothreshold[[x]], 
                genomicpercentidentity = genomicpercentidentity[[x]])
            if (!as.logical(keepmultihits[[x]])) {
                psl.x <- psl.x[!psl.x$isMultiHit, ]
            }
            psl.x <- clusterSites(psl.rd = psl.x, windowSize = clustersiteswithin[[x]])
            if (as.logical(keepmultihits[[x]])) {
                psl.x <- otuSites(psl.rd = psl.x)
            }
            psl.x
        }
    names(psl.hits) <- names(psl)
    message("Adding PSL hits back to the object.")
    sampleInfo <- addFeature(sampleInfo, sector = NULL, samplename = names(psl.hits), 
        feature = "psl", value = psl.hits)
    message("Adding sites back to the object.")
    psl.hits <- sapply(psl.hits, function(x) x[with(x, clusterTopHit & 
        pass.allQC & !isMultiHit), ])
    sampleInfo <- addFeature(sampleInfo, sector = NULL, samplename = names(psl.hits), 
        feature = "sites", value = psl.hits)
    cleanit <- gc()
    sampleInfo$callHistory <- append(sampleInfo$callHistory, 
        match.call())
    return(sampleInfo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
