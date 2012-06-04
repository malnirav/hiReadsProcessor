\name{summary.simple}
\alias{summary.simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.simple(sampleInfo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfo}{
%%     ~~Describe \code{sampleInfo} here~~
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
function (sampleInfo) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    cat("Total sectors:", paste(names(sampleInfo$sectors), collapse = ","), 
        "\n")
    for (sector in names(sampleInfo$sectors)) {
        cat(rep("~", 40), "\n")
        cat("Sector:", sector, "\n")
        for (metaD in c("decoded", "primed", "LTRed", "linkered", 
            "psl", "sites")) {
            res <- sapply(extractFeature(sampleInfo, sector = sector, 
                feature = metaD)[[sector]], length)
            if (length(res) > 0) {
                cat("\tTotal", metaD, ":", sum(res), "\n")
                cat("\t\tBy Sample:", paste(names(res), res, 
                  sep = ":", collapse = "\n\t\t\t   "), "\n")
            }
            else {
                cat("\tNo", metaD, "attribute found....skipping\n")
            }
            cat("\t", rep("-", 25), "\n")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
