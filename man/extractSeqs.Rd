\name{extractSeqs}
\alias{extractSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractSeqs(sampleInfo, sector = NULL, samplename = NULL, feature = "genomic", trim = TRUE, minReadLength = 1, sideReturn = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfo}{
%%     ~~Describe \code{sampleInfo} here~~
}
  \item{sector}{
%%     ~~Describe \code{sector} here~~
}
  \item{samplename}{
%%     ~~Describe \code{samplename} here~~
}
  \item{feature}{
%%     ~~Describe \code{feature} here~~
}
  \item{trim}{
%%     ~~Describe \code{trim} here~~
}
  \item{minReadLength}{
%%     ~~Describe \code{minReadLength} here~~
}
  \item{sideReturn}{
%%     ~~Describe \code{sideReturn} here~~
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
function (sampleInfo, sector = NULL, samplename = NULL, feature = "genomic", 
    trim = TRUE, minReadLength = 1, sideReturn = NULL) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    if (is.null(feature)) {
        stop("Please define a feature to extract.")
    }
    res <- getSectorsForSamples(sampleInfo, sector, samplename)
    sectors <- res[["sectors"]]
    samplenames <- res[["samplenames"]]
    if (feature == "unDecoded") {
        sapply(sectors, function(y) {
            allmetadata <- metadata(sampleInfo$sectors[[y]])
            if ("unDecodedSeqs" \%in\% names(allmetadata)) {
                allmetadata$unDecodedSeqs
            }
            else {
                message("No unDecoded attribute found for the supplied sampleInfo object and sector ", 
                  y)
            }
        })
    }
    else {
        res <- sapply(sectors, function(y) {
            sapply(samplenames[[y]], function(x, y) {
                decoded <- sampleInfo$sectors[[y]]$samples[[x]]$decoded
                if (feature != "decoded") {
                  assign(gsub("!", "", feature), sampleInfo$sectors[[y]]$samples[[x]][[gsub("!", 
                    "", feature)]])
                }
                if (feature == "decoded") {
                  decoded
                }
                else if (feature \%in\% c("genomic", "genomicLinkered")) {
                  primed <- sampleInfo$sectors[[y]]$samples[[x]]$primed
                  LTRed <- sampleInfo$sectors[[y]]$samples[[x]]$LTRed
                  linkered <- sampleInfo$sectors[[y]]$samples[[x]]$linkered
                  if (is.null(LTRed)) {
                    warning("LTRed information not found...using primer end as starting boundary.", 
                      immediate. = TRUE)
                  }
                  if (any(is.null(primed), is.null(linkered))) {
                    message("No sequences found for requested feature (", 
                      feature, ") for sample: ", x, "...skipping.")
                  }
                  else {
                    if (trim) {
                      if (is.null(LTRed)) {
                        LTRed <- primed
                      }
                      starts <- end(LTRed)
                      starts.name <- names(LTRed)
                      ends <- width(decoded)
                      ends.name <- names(decoded)
                      stopifnot(identical(starts.name, ends.name[ends.name \%in\% 
                        starts.name]))
                      coords <- IRanges(start = starts + 1, end = ends[ends.name \%in\% 
                        starts.name], names = starts.name)
                      ends <- start(linkered)
                      ends.name <- names(linkered)
                      end(coords[ends.name]) <- ends - 1
                      if (feature == "genomicLinkered") {
                        coords <- coords[names(coords) \%in\% names(linkered)]
                      }
                      seqs <- trimSeqs(decoded, coords)
                      seqs[width(seqs) >= minReadLength]
                    }
                    else {
                      if (feature == "genomicLinkered") {
                        decoded[names(decoded) \%in\% names(linkered)]
                      }
                      else {
                        if (is.null(LTRed)) {
                          decoded[names(decoded) \%in\% names(primed)]
                        }
                        else {
                          decoded[names(decoded) \%in\% names(LTRed)]
                        }
                      }
                    }
                  }
                }
                else {
                  notFeature <- grepl("!", feature)
                  if (notFeature) {
                    toreturn <- decoded[!names(decoded) \%in\% 
                      names(get(gsub("!", "", feature)))]
                    if (length(toreturn)) {
                      toreturn
                    }
                    else {
                      message("No sequences found for requested feature (", 
                        feature, ") for sample: ", x, "...skipping.")
                    }
                  }
                  else {
                    if (is.null(get(feature))) {
                      message("No sequences found for requested feature (", 
                        feature, ") for sample: ", x, "...skipping.")
                    }
                    else {
                      res.seq <- decoded[names(decoded) \%in\% 
                        names(get(feature))]
                      if (trim) {
                        if (is.null(sideReturn)) {
                          sidetype <- ifelse(grepl("primerID", 
                            feature, ignore.case = TRUE), "middle", 
                            ifelse(grepl("primed|LTRed", feature, 
                              ignore.case = TRUE), "left", ifelse(grepl("linkered", 
                              feature, ignore.case = TRUE), "right", 
                              "middle")))
                        }
                        else {
                          sidetype <- tolower(sideReturn)
                        }
                        offByLength <- ifelse(sidetype == "middle", 
                          0, 1)
                        seqs <- trimSeqs(res.seq, get(feature), 
                          side = sidetype, offBy = offByLength)
                        seqs[width(seqs) >= minReadLength]
                      }
                      else {
                        res.seq
                      }
                    }
                  }
                }
            }, y = y)
        }, simplify = FALSE)
        lengthTest <- lapply(lapply(lapply(res, function(x) sapply(x, 
            length)), ">", 0), which)
        mapply(function(x, y) x[y], res, lengthTest, SIMPLIFY = FALSE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
