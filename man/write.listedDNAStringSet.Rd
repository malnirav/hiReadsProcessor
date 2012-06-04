\name{write.listedDNAStringSet}
\alias{write.listedDNAStringSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.listedDNAStringSet(dnaSet, filePath = ".", filePrefix = "processed", prependSamplenames = TRUE, parallel = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dnaSet}{
%%     ~~Describe \code{dnaSet} here~~
}
  \item{filePath}{
%%     ~~Describe \code{filePath} here~~
}
  \item{filePrefix}{
%%     ~~Describe \code{filePrefix} here~~
}
  \item{prependSamplenames}{
%%     ~~Describe \code{prependSamplenames} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
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
function (dnaSet, filePath = ".", filePrefix = "processed", prependSamplenames = TRUE, 
    parallel = TRUE) 
{
    stopifnot(class(dnaSet) == "list")
    if (!parallel) {
        registerDoSEQ()
    }
    for (sector in names(dnaSet)) {
        foreach(i = iter(names(dnaSet[[sector]]), .inorder = FALSE, 
            .errorhandling = "stop", .packages = "Biostrings", 
            .export = c("dnaSet", "sector", "filePath", "filePrefix", 
                "prependSamplenames"))) \%dopar\% {
            outputSeqs <- dnaSet[[sector]][[i]]
            if (length(outputSeqs) > 0) {
                if (is.null(names(outputSeqs))) {
                  message("No names attribute found for ", i, 
                    " ... using artifically generated names")
                  names(outputSeqs) <- paste("read", 1:length(outputSeqs), 
                    sep = "-")
                }
                filename <- paste(filePath, paste(filePrefix, 
                  i, "fa", sep = "."), sep = "/")
                if (prependSamplenames) {
                  names(outputSeqs) <- paste(i, names(outputSeqs), 
                    sep = "-")
                }
                write.XStringSet(outputSeqs, file = filename, 
                  format = "fasta", append = TRUE, width = 20000)
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
