\name{addFeature}
\alias{addFeature}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
addFeature(sampleInfo, sector = NULL, samplename = NULL, feature = NULL, value = NULL)
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
  \item{value}{
%%     ~~Describe \code{value} here~~
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
function (sampleInfo, sector = NULL, samplename = NULL, feature = NULL, 
    value = NULL) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    if (is.null(feature)) {
        stop("Please define the new feature(s) to add.")
    }
    if (is.null(value)) {
        stop("Please define the value parameter.")
    }
    if (is.null(samplename)) {
        stop("Please define the samplename to add the features to.")
    }
    if (!all(samplename \%in\% names(value))) {
        stop("Not all samplename(s) are found in value parameter")
    }
    res <- getSectorsForSamples(sampleInfo, sector, samplename)
    sectors <- res[["sectors"]]
    samplenames <- res[["samplenames"]]
    for (y in sectors) {
        for (x in samplenames[[y]]) {
            cat(".")
            sampleInfo$sectors[[y]]$samples[[x]][[feature]] <- value[[x]]
        }
    }
    return(sampleInfo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
