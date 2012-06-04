\name{splitSeqsToFiles}
\alias{splitSeqsToFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splitSeqsToFiles(x, totalFiles = 4, suffix = "tempy", filename = "queryFile.fa")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{totalFiles}{
%%     ~~Describe \code{totalFiles} here~~
}
  \item{suffix}{
%%     ~~Describe \code{suffix} here~~
}
  \item{filename}{
%%     ~~Describe \code{filename} here~~
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
function (x, totalFiles = 4, suffix = "tempy", filename = "queryFile.fa") 
{
    if (is.atomic(x)) {
        message("Splitting file ", x)
        totalSeqs <- length(fasta.info(x, use.names = FALSE))
        chunks <- round(totalSeqs/totalFiles)
        chunks <- ifelse(chunks > 0, chunks, totalSeqs)
        starts <- seq(0, totalSeqs, by = chunks)
        for (skippy in starts[starts != totalSeqs]) {
            filename.out <- paste(x, skippy, suffix, sep = ".")
            query.tmp <- read.BStringSet(x, nrec = chunks, skip = skippy)
            write.XStringSet(query.tmp, file = filename.out, 
                format = "fasta")
        }
        return(list.files(path = dirname(x), pattern = paste(basename(x), 
            ".*", suffix, "$", sep = ""), full.names = TRUE))
    }
    else if (class(x) == "DNAStringSet") {
        message("Splitting Reads.")
        totalSeqs <- length(x)
        chunks <- round(totalSeqs/totalFiles)
        starts <- seq(1, totalSeqs, by = chunks)
        stops <- unique(c(seq(chunks, totalSeqs, by = chunks), 
            totalSeqs))
        stopifnot(length(starts) == length(stops))
        for (skippy in 1:length(starts)) {
            filename.out <- paste(filename, skippy, suffix, sep = ".")
            write.XStringSet(x[starts[skippy]:stops[skippy]], 
                file = filename.out, format = "fasta")
        }
        return(list.files(path = ".", pattern = paste(filename, 
            ".*", suffix, "$", sep = ""), full.names = TRUE))
    }
    else {
        stop("Dont know what is supplied in parameter x.")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
