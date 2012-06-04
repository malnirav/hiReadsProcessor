\name{splitByBarcode}
\alias{splitByBarcode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splitByBarcode(barcodesSample, dnaSet, trimFrom = NULL, showStats = FALSE, returnUnmatched = FALSE, dereplicate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{barcodesSample}{
%%     ~~Describe \code{barcodesSample} here~~
}
  \item{dnaSet}{
%%     ~~Describe \code{dnaSet} here~~
}
  \item{trimFrom}{
%%     ~~Describe \code{trimFrom} here~~
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
function (barcodesSample, dnaSet, trimFrom = NULL, showStats = FALSE, 
    returnUnmatched = FALSE, dereplicate = FALSE) 
{
    if (is.null(barcodesSample) | length(barcodesSample) == 0) {
        stop("No barcodes to samples association vector provided in parameter barcodesSample.")
    }
    if (is.null(dnaSet) | length(dnaSet) == 0) {
        stop("No sequences provided in parameter dnaSet.")
    }
    if (is.null(names(dnaSet))) {
        stop("No names attribute found in dnaSet object")
    }
    message("Using following schema for barcode to sample associations")
    print(as.data.frame(barcodesSample))
    barcodelen <- unique(nchar(names(barcodesSample)))
    seqbarcodes <- substr(dnaSet, 1, barcodelen)
    good.row <- seqbarcodes \%in\% names(barcodesSample)
    if (!any(good.row)) {
        stop("No matching barcoded sequences found on this quadrant.")
    }
    sampleNames <- barcodesSample[seqbarcodes[good.row]]
    deflines <- sub("^(\\S+) .+$", "\\1", names(dnaSet)[good.row], 
        perl = T)
    if (is.null(trimFrom)) {
        trimFrom <- barcodelen + 1
    }
    unmatched <- DNAStringSet(dnaSet[!good.row])
    dnaSet <- DNAStringSet(dnaSet[good.row], start = trimFrom)
    names(dnaSet) <- deflines
    if (showStats) {
        message("Number of Sequences with no matching barcode: ", 
            as.numeric(table(good.row)["FALSE"]))
        message("Number of Sequences decoded:")
        print(as.data.frame(table(sampleNames)))
    }
    dnaSet <- split(dnaSet, as.character(sampleNames))
    if (dereplicate) {
        message("Dereplicating reads.")
        dnaSet <- lapply(dnaSet, dereplicateReads)
    }
    if (returnUnmatched) {
        dnaSet <- c(dnaSet, unDecodedSeqs = unmatched)
    }
    return(dnaSet)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
