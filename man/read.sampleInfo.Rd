\name{read.sampleInfo}
\alias{read.sampleInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.sampleInfo(sampleInfoPath = NULL, splitBySector = TRUE, interactive = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfoPath}{
%%     ~~Describe \code{sampleInfoPath} here~~
}
  \item{splitBySector}{
%%     ~~Describe \code{splitBySector} here~~
}
  \item{interactive}{
%%     ~~Describe \code{interactive} here~~
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
function (sampleInfoPath = NULL, splitBySector = TRUE, interactive = TRUE) 
{
    if (is.null(sampleInfoPath)) {
        stop("No sample information file path provided.")
    }
    sampleInfoPath <- normalizePath(sampleInfoPath, mustWork = TRUE)
    requiredCols <- c("sector", "barcode", "primerltrsequence", 
        "samplename", "sampledescription", "gender", "species", 
        "freeze", "linkersequence", "restrictionenzyme")
    metaDataCols <- c(ltrbitsequence = "", ltrbitidentity = 1, 
        primerltridentity = 0.85, linkeridentity = 0.55, primeridinlinker = FALSE, 
        primeridinlinkeridentity1 = 0.75, primeridinlinkeridentity2 = 0.5, 
        celltype = "", user = Sys.getenv("USER"), startwithin = 3, 
        alignratiothreshold = 0.7, clustersiteswithin = 5, keepmultihits = TRUE, 
        genomicpercentidentity = 0.98, processingdate = format(Sys.time(), 
            "\%Y-\%m-\%d "))
    if (grepl(".xls.?$", sampleInfoPath)) {
        sampleInfo <- unique(read.xlsx(sampleInfoPath, sheetIndex = 1, 
            stringsAsFactors = FALSE))
    }
    else {
        sampleInfo <- unique(read.delim(sampleInfoPath, stringsAsFactors = FALSE))
    }
    names(sampleInfo) <- tolower(gsub("\\.|-|_", "", names(sampleInfo)))
    ColsNotThere <- !requiredCols \%in\% names(sampleInfo)
    if (any(ColsNotThere)) {
        absentCols <- requiredCols[ColsNotThere]
        info <- paste("Following required column(s) is absent from the Sample Info file:", 
            paste(absentCols, sep = "", collapse = ", "), sep = "")
        stop(info)
    }
    metaColsNotThere <- !names(metaDataCols) \%in\% names(sampleInfo)
    if (any(metaColsNotThere)) {
        sampleInfo <- cbind(sampleInfo, as.data.frame(t(metaDataCols[metaColsNotThere]), 
            stringsAsFactors = FALSE))
    }
    for (column in c("sector", "barcode", "primerltrsequence", 
        "ltrbitsequence", "samplename", "linkersequence", "restrictionenzyme")) {
        sampleInfo[, column] <- gsub(" ", "", sampleInfo[, column])
        if (column \%in\% c("barcode", "primerltrsequence", "ltrbitsequence", 
            "linkersequence", "restrictionenzyme")) {
            sampleInfo[, column] <- toupper(sampleInfo[, column])
        }
    }
    ltrbitTest <- nchar(sampleInfo$ltrbitsequence) == 0 | sampleInfo$ltrbitsequence == 
        ""
    if (any(ltrbitTest)) {
        tofix <- which(ltrbitTest)
        if (interactive) {
            message("LTR bit not found for ", length(tofix), 
                " samples. Use last 7 bases of the LTR primer as the LTR bit? (y/n)")
            choice <- scan(what = character(0), n = 1, quiet = TRUE, 
                multi.line = FALSE)
        }
        else {
            message("LTR bit not found for ", length(tofix), 
                " samples. Using last 7 bases of the LTR primer as the LTR bit.")
            choice <- "y"
        }
        if (tolower(choice) == "y") {
            sampleInfo$ltrbitsequence <- substr(sampleInfo$primerltrsequence, 
                nchar(sampleInfo$primerltrsequence) - 6, nchar(sampleInfo$primerltrsequence))
            sampleInfo$primerltrsequence <- substr(sampleInfo$primerltrsequence, 
                1, nchar(sampleInfo$primerltrsequence) - 7)
        }
        else {
            warning("No LTR bit sequence found for following samples: ", 
                paste(sampleInfo$samplename[tofix], sep = "", 
                  collapse = ", "), immediate. = TRUE)
        }
    }
    samplenametest <- nchar(sampleInfo$samplename) == 0 | sampleInfo$samplename == 
        ""
    if (any(samplenametest)) {
        stop("No sample names found in following rows of the sample information file ", 
            sampleInfoPath, " : ", paste(which(samplenametest), 
                sep = "", collapse = ", "))
    }
    sectortest <- nchar(sampleInfo$sector) == 0 | sampleInfo$sector == 
        "" | is.na(sampleInfo$sector)
    if (any(sectortest)) {
        tofix <- which(sectortest)
        if (interactive) {
            message("Sector information not found for ", length(tofix), 
                " samples. Which sector are they from? (1,2,4,etc)")
            choice <- scan(what = character(0), quiet = TRUE, 
                multi.line = FALSE)
        }
        else {
            message("Sector information not found for ", length(tofix), 
                " samples. Assuming they are from sector 1.")
            choice <- "1"
        }
        if (length(choice) > 0) {
            sampleInfo$sector[tofix] <- unlist(strsplit(choice, 
                ","))
        }
        else {
            stop("No Sector information found for following samples: ", 
                paste(sampleInfo$samplename[tofix], sep = "", 
                  collapse = ", "))
        }
    }
    sampleInfo$sector <- gsub("\\.0$", "", as.character(sampleInfo$sector))
    sampleSectorTest <- table(paste(sampleInfo$samplename, sampleInfo$sector))
    if (any(sampleSectorTest > 1)) {
        stop("Duplicate sample names found on same quadrant in the sample information file ", 
            sampleInfoPath, " : ", paste(sampleSectorTest[sampleSectorTest > 
                1], sep = "", collapse = ", "))
    }
    if (splitBySector) {
        sampleInfo <- SimpleList(split(sampleInfo, sampleInfo$sector))
        for (sector in 1:length(sampleInfo)) {
            sampleData <- do.call(SimpleList, split(sampleInfo[[sector]], 
                sampleInfo[[sector]]$samplename))
            for (sample.i in 1:length(sampleData)) {
                sampleData[[sample.i]] <- do.call(SimpleList, 
                  as.list(sampleData[[sample.i]]))
            }
            sampleInfo[[sector]] <- SimpleList(samples = sampleData)
        }
    }
    return(sampleInfo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
