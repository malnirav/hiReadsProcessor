\name{removeReadsWithNs}
\alias{removeReadsWithNs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
removeReadsWithNs(dnaSet, maxNs = 5, consecutive = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dnaSet}{
%%     ~~Describe \code{dnaSet} here~~
}
  \item{maxNs}{
%%     ~~Describe \code{maxNs} here~~
}
  \item{consecutive}{
%%     ~~Describe \code{consecutive} here~~
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
function (dnaSet, maxNs = 5, consecutive = TRUE) 
{
    if (consecutive) {
        good.row <- grepl(paste(rep("N", maxNs), collapse = ""), 
            dnaSet, fixed = TRUE)
    }
    else {
        res <- alphabetFrequency(dnaSet)
        good.row <- res[, "N"] <= maxNs
    }
    return(dnaSet[good.row])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
