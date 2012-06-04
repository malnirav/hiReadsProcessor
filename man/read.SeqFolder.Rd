\name{read.SeqFolder}
\alias{read.SeqFolder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.SeqFolder(SequencingFolderPath = NULL, sampleInfoFilePath = NULL, seqfilePattern = "\\.TCA.454Reads.fna$", interactive = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SequencingFolderPath}{
%%     ~~Describe \code{SequencingFolderPath} here~~
}
  \item{sampleInfoFilePath}{
%%     ~~Describe \code{sampleInfoFilePath} here~~
}
  \item{seqfilePattern}{
%%     ~~Describe \code{seqfilePattern} here~~
}
  \item{interactive}{
%%     ~~Describe \code{interactive} here~~
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
function (SequencingFolderPath = NULL, sampleInfoFilePath = NULL, 
    seqfilePattern = "\\.TCA.454Reads.fna$", interactive = TRUE) 
{
    if (is.null(SequencingFolderPath)) {
        stop("No Sequencing Folder Path provided.")
    }
    SequencingFolderPath <- normalizePath(SequencingFolderPath, 
        mustWork = TRUE)
    seqFilePaths <- list.files(path = SequencingFolderPath, recursive = TRUE, 
        full.names = TRUE, pattern = seqfilePattern)
    if (length(seqFilePaths) == 0) {
        stop(paste("No files found in the folder", SequencingFolderPath, 
            "matching following pattern:", seqfilePattern))
    }
    if (is.null(sampleInfoFilePath)) {
        possibleFiles <- list.files(path = SequencingFolderPath, 
            recursive = TRUE, full.names = TRUE, pattern = ".*sampleinfo.+", 
            ignore.case = TRUE)
        if (length(possibleFiles) == 0) {
            stop(paste("No sample information file found in folder :", 
                SequencingFolderPath))
        }
        else {
            if (interactive & length(possibleFiles) > 1) {
                message("Please choose a sample information file to read the meta data from:\n", 
                  paste(1:length(possibleFiles), possibleFiles, 
                    sep = ": ", collapse = "\n"))
                choice <- scan(what = integer(0), n = 1, quiet = TRUE, 
                  multi.line = FALSE)
            }
            else {
                choice <- 1
            }
            message("Choosing ", possibleFiles[choice], " as sample information file.")
        }
        sampleInfoFilePath <- possibleFiles[choice]
    }
    sampleInfo <- read.sampleInfo(sampleInfoFilePath, interactive = interactive)
    if (length(sampleInfo) != length(seqFilePaths)) {
        warning("Number of sectors (", length(sampleInfo), ") in sample information file does not match number of sector files (", 
            length(seqFilePaths), ") found in the folder.")
    }
    return(SimpleList(SequencingFolderPath = SequencingFolderPath, 
        seqFilePaths = seqFilePaths, seqfilePattern = seqfilePattern, 
        sampleInfoFilePath = sampleInfoFilePath, sectors = sampleInfo, 
        callHistory = match.call()))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
