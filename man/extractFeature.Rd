\name{extractFeature}
\alias{extractFeature}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractFeature(sampleInfo, sector = NULL, samplename = NULL, feature = NULL)
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
function (sampleInfo, sector = NULL, samplename = NULL, feature = NULL) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    if (is.null(feature)) {
        stop("Please define a feature to extract.")
    }
    res <- getSectorsForSamples(sampleInfo, sector, samplename)
    sectors <- res[["sectors"]]
    samplenames <- res[["samplenames"]]
    res <- sapply(sectors, function(y) {
        sapply(samplenames[[y]], function(x, y) {
            if (feature == "metadata") {
                paste(names(sampleInfo$sectors[[y]]$samples[[x]]), 
                  collapse = ", ")
            }
            else {
                res <- sampleInfo$sectors[[y]]$samples[[x]][[feature]]
                if (class(res) == "factor") {
                  if (!any(is.na(suppressWarnings(as.numeric(levels(res)))))) {
                    as.numeric(as.character(res))
                  }
                  else {
                    as.character(res)
                  }
                }
                else if (class(res) == "character") {
                  if (!any(is.na(suppressWarnings(as.numeric(res))))) {
                    as.numeric(res)
                  }
                  else {
                    res
                  }
                }
                else {
                  res
                }
            }
        }, y = y)
    }, simplify = FALSE)
    lengthTest <- lapply(lapply(lapply(res, function(x) sapply(x, 
        length)), ">", 0), which)
    mapply(function(x, y) x[y], res, lengthTest, SIMPLIFY = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
