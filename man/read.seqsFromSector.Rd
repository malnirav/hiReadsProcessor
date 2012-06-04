\name{read.seqsFromSector}
\alias{read.seqsFromSector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.seqsFromSector(seqFilePath = NULL, sector = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqFilePath}{
%%     ~~Describe \code{seqFilePath} here~~
}
  \item{sector}{
%%     ~~Describe \code{sector} here~~
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
function (seqFilePath = NULL, sector = 1) 
{
    if (is.null(seqFilePath)) {
        stop("Missing seqFilePath!")
    }
    if (class(seqFilePath) == "SimpleList") {
        filePath <- normalizePath(grep(paste(sector, seqFilePath$seqfilePattern, 
            sep = ""), seqFilePath$seqFilePaths, value = TRUE), 
            mustWork = TRUE)
        if (!length(filePath) > 0) {
            stop("No sequence file found for sector: ", sector, 
                " in seqFilePath variable (", seqFilePath$seqFilePaths, 
                ") using pattern (", paste(sector, seqFilePath$seqfilePattern, 
                  sep = ""), ")")
        }
        seqFilePath <- filePath
    }
    message("Reading ", seqFilePath)
    dnaSet <- read.DNAStringSet(normalizePath(seqFilePath, mustWork = TRUE), 
        format = ifelse(grepl("fastq$", seqFilePath), "fastq", 
            "fasta"))
    if (length(dnaSet) == 0) {
        stop("No sequences found")
    }
    return(dnaSet)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
