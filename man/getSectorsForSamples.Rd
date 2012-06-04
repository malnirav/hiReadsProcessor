\name{getSectorsForSamples}
\alias{getSectorsForSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSectorsForSamples(sampleInfo, sector = NULL, samplename = NULL, returnDf = FALSE)
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
  \item{returnDf}{
%%     ~~Describe \code{returnDf} here~~
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
function (sampleInfo, sector = NULL, samplename = NULL, returnDf = FALSE) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    if (is.null(sector)) {
        sectors <- names(sampleInfo$sectors)
    }
    else {
        sectors <- sector
    }
    samplenames <- sapply(sectors, function(x) names(sampleInfo$sectors[[x]]$samples), 
        simplify = FALSE)
    if (is.null(samplename)) {
        samplename <- unlist(samplenames)
    }
    if (!all(samplename \%in\% unlist(samplenames))) {
        stop("Following sample(s) do not exist on given sector(s) (", 
            paste(sectors, collapse = ", "), ") in the supplied sampleInfo object: ", 
            paste(samplename[!samplename \%in\% unlist(samplenames)], 
                collapse = ", "))
    }
    sectors <- names(which(unlist(lapply(lapply(samplenames, 
        "\%in\%", samplename), any))))
    if (returnDf) {
        return(do.call(rbind, lapply(sectors, function(x) {
            data.frame(samplename = samplenames[[x]][samplenames[[x]] \%in\% 
                samplename], sector = x, stringsAsFactors = FALSE)
        })))
    }
    else {
        samplenames <- sapply(sectors, function(x) samplenames[[x]][samplenames[[x]] \%in\% 
            samplename], simplify = FALSE)
        return(list(sectors = sectors, samplenames = samplenames))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
