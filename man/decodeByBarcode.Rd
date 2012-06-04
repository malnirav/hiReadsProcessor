\name{decodeByBarcode}
\alias{decodeByBarcode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
decodeByBarcode(sampleInfo, sector = NULL, dnaSet = NULL, showStats = FALSE, returnUnmatched = FALSE, dereplicate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfo}{
%%     ~~Describe \code{sampleInfo} here~~
}
  \item{sector}{
%%     ~~Describe \code{sector} here~~
}
  \item{dnaSet}{
%%     ~~Describe \code{dnaSet} here~~
}
  \item{showStats}{
%%     ~~Describe \code{showStats} here~~
}
  \item{returnUnmatched}{
%%     ~~Describe \code{returnUnmatched} here~~
}
  \item{dereplicate}{
%%     ~~Describe \code{dereplicate} here~~
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
function (sampleInfo, sector = NULL, dnaSet = NULL, showStats = FALSE, 
    returnUnmatched = FALSE, dereplicate = FALSE) 
{
    if (class(sampleInfo) == "SimpleList") {
        if (is.null(sector)) {
            stop("No sector provided in parameter sector.")
        }
        sectors <- sector <- as.character(sector)
        if (length(sectors) == 1 & tolower(sectors) == "all") {
            sectors <- names(sampleInfo$sectors)
        }
        if (any(!sectors \%in\% names(sampleInfo$sectors))) {
            stop("Following sectors not found in names(sampleInfo$sectors): ", 
                sectors[!sectors \%in\% names(sampleInfo$sectors)])
        }
        for (sector in sectors) {
            message("Decoding sector: ", sector)
            sampleBarcodes <- toupper(extractFeature(sampleInfo, 
                sector = sector, feature = "barcode")[[sector]])
            barcodesSample <- c()
            barcodesSample[as.character(sampleBarcodes)] <- names(sampleBarcodes)
            if (length(table(nchar(as.character(sampleBarcodes)))) > 
                1) {
                stop("Multiple barcode lengths found.")
            }
            realbarcodelen <- unique(nchar(as.character(sampleBarcodes)))
            if (any(table(as.character(sampleBarcodes)) > 1)) {
                message("Duplicate barcode found on this sector.\nPlease choose from one of the options below\n 1: Pick first few bases of primer for tiebreaker? (This could be dangerous if the sequencing run has too many errors!)\n 2: Use the last sample associated with the duplicate as the primary sample?\n 3: Do not do anything.")
                choice <- scan(what = integer(0), n = 1, quiet = TRUE, 
                  multi.line = FALSE)
                if (choice == 1) {
                  message("Enter # of bases to use from primer:")
                  howmany <- scan(what = integer(0), n = 1, quiet = TRUE, 
                    multi.line = FALSE)
                  samplePrimers <- toupper(extractFeature(sampleInfo, 
                    sector = sector, feature = "primerltrsequence")[[sector]])
                  newBarcodes <- toupper(paste(sampleBarcodes, 
                    substr(samplePrimers, 1, howmany), sep = ""))
                  if (any(table(newBarcodes) > 1)) {
                    stop("Tiebreaking failed...try choosing high number of bases from primer??? \n", 
                      paste(names(which(table(newBarcodes) > 
                        1)), collapse = ", "))
                  }
                  barcodesSample <- c()
                  barcodesSample[newBarcodes] <- names(sampleBarcodes)
                }
                else if (choice == 2) {
                  message("Overwriting duplicate samples associated with the same barcode...")
                }
                else {
                  stop("Aborting due to duplicate barcode found on this sector")
                }
            }
            dnaSet <- read.seqsFromSector(sampleInfo, sector)
            dnaSet <- splitByBarcode(barcodesSample, dnaSet, 
                trimFrom = realbarcodelen + 1, showStats = showStats, 
                returnUnmatched = returnUnmatched, dereplicate = dereplicate)
            for (samplename in names(dnaSet)) {
                if (samplename == "unDecodedSeqs") {
                  metadata(sampleInfo$sectors[[sector]]) <- append(metadata(sampleInfo$sectors[[sector]]), 
                    list(unDecodedSeqs = dnaSet[[samplename]]))
                }
                else {
                  sampleInfo$sectors[[sector]]$samples[[samplename]]$decoded <- dnaSet[[samplename]]
                }
            }
            metadata(sampleInfo$sectors[[sector]]) <- append(metadata(sampleInfo$sectors[[sector]]), 
                list(decodedBy = barcodesSample))
        }
        sampleInfo$callHistory <- append(sampleInfo$callHistory, 
            match.call())
        decoded <- sampleInfo
        cleanit <- gc()
    }
    else {
        decoded <- splitByBarcode(sampleInfo, dnaSet, trimFrom = NULL, 
            showStats = showStats, returnUnmatched = returnUnmatched, 
            dereplicate = dereplicate)
        cleanit <- gc()
    }
    return(decoded)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
