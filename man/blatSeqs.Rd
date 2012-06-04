\name{blatSeqs}
\alias{blatSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
blatSeqs(query = NULL, subject = NULL, standaloneBlat = TRUE, port = NULL, host = "localhost", parallel = TRUE, gzipResults = TRUE, blatParameters = c(minIdentity = 70, minScore = 5, stepSize = 8, tileSize = 10, repMatch = 112312, dots = 50, q = "dna", t = "dna"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{
%%     ~~Describe \code{query} here~~
}
  \item{subject}{
%%     ~~Describe \code{subject} here~~
}
  \item{standaloneBlat}{
%%     ~~Describe \code{standaloneBlat} here~~
}
  \item{port}{
%%     ~~Describe \code{port} here~~
}
  \item{host}{
%%     ~~Describe \code{host} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{gzipResults}{
%%     ~~Describe \code{gzipResults} here~~
}
  \item{blatParameters}{
%%     ~~Describe \code{blatParameters} here~~
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
function (query = NULL, subject = NULL, standaloneBlat = TRUE, 
    port = NULL, host = "localhost", parallel = TRUE, gzipResults = TRUE, 
    blatParameters = c(minIdentity = 70, minScore = 5, stepSize = 8, 
        tileSize = 10, repMatch = 112312, dots = 50, q = "dna", 
        t = "dna")) 
{
    suppressWarnings(blatOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", 
        "\\1", grep("\\s+-.+=", system("blat", intern = TRUE), 
            value = TRUE))))
    suppressWarnings(gfClientOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", 
        "\\1", grep("\\s+-.+=", system("gfclient", intern = TRUE), 
            value = TRUE))))
    if (!standaloneBlat) {
        message("Using gfClient protocol to perform BLAT.")
        if (is.null(port)) {
            stop("The port paramter is empty. Please define the port used to start gfServer with")
        }
    }
    if (is.null(subject)) {
        stop("The subject parameter is empty. Please supply subject sequences or a path to 2bit or nib files to serve as reference/target")
    }
    else {
        subjectFile <- NULL
        if (is.atomic(subject)) {
            if (any(grepl("\\.2bit$|\\.nib$", subject))) {
                if (standaloneBlat) {
                  stop("Standalone BLAT cannot be used when subject is an indexed nib or 2bit file.")
                }
                indexFileDir <- dirname(subject)
                subjectFile <- list.files(path = indexFileDir, 
                  pattern = basename(subject), full.names = TRUE)
            }
            else {
                subject <- DNAStringSet(subject)
            }
        }
        if (is.null(subjectFile)) {
            if (is.null(names(subject))) {
                names(subject) <- paste("subject", 1:length(subject))
            }
            filename.seq <- "subjectFile.fa.tempyS"
            write.XStringSet(subject, file = filename.seq, format = "fasta")
            subjectFile <- filename.seq
        }
    }
    if (is.null(query)) {
        stop("The query parameter is empty. Please supply reads to be aligned")
    }
    else {
        queryFiles <- NULL
        if (is.atomic(query)) {
            if (any(grepl("\\.fna$|\\.fa$|\\*", query))) {
                queryFiles <- list.files(path = dirname(query), 
                  pattern = basename(query), full.names = TRUE)
                if (parallel) {
                  queryFiles <- unlist(sapply(queryFiles, function(f) splitSeqsToFiles(f, 
                    getDoParWorkers(), "tempyQ")), use.names = FALSE)
                }
            }
            else {
                query <- DNAStringSet(query)
            }
        }
        if (is.null(queryFiles)) {
            if (is.null(names(query))) {
                names(query) <- paste("read", 1:length(query))
            }
            if (parallel) {
                queryFiles <- splitSeqsToFiles(query, getDoParWorkers(), 
                  "tempyQ")
            }
            else {
                queryFiles <- "queryFile.fa.tempyQ"
                write.XStringSet(query, file = queryFiles, format = "fasta")
            }
        }
    }
    if (!parallel) {
        registerDoSEQ()
    }
    queryFiles <- as.character(queryFiles)
    subjectFile <- as.character(subjectFile)
    if (standaloneBlat) {
        blatOpts <- blatParameters[names(blatParameters) \%in\% 
            blatOpts]
        stopifnot(length(subjectFile) == 1)
        filenames <- foreach(x = iter(queryFiles), .inorder = FALSE, 
            .export = c("blatOpts", "subjectFile", "gzipResults")) \%dopar\% 
            {
                filename.psl <- paste(x, "psl", sep = ".")
                cmd <- paste("blat", paste(paste("-", names(blatOpts), 
                  sep = ""), blatOpts, collapse = " ", sep = "="), 
                  "-noHead", subjectFile, x, filename.psl)
                message(cmd)
                system(cmd)
                if (grepl("\\.tempyQ$", x)) {
                  system(paste("rm", x))
                }
                if (gzipResults) {
                  system(paste("gzip", filename.psl))
                  filename.psl <- paste(filename.psl, "gz", sep = ".")
                }
                filename.psl
            }
        if (grepl("\\.tempyS$", subjectFile)) {
            system(paste("rm", subjectFile))
        }
    }
    else {
        gfClientOpts <- blatParameters[names(blatParameters) \%in\% 
            gfClientOpts]
        stopifnot(length(subjectFile) > 0)
        filenames <- foreach(x = iter(queryFiles), .inorder = FALSE, 
            .export = c("gfClientOpts", "host", "port", "indexFileDir", 
                "gzipResults")) \%dopar\% {
            filename.psl <- paste(x, "psl", sep = ".")
            cmd <- paste("gfClient", paste(paste("-", names(gfClientOpts), 
                sep = ""), gfClientOpts, collapse = " ", sep = "="), 
                "-nohead", host, port, "/", x, filename.psl)
            message(cmd)
            system(cmd)
            if (grepl("\\.tempyQ$", x)) {
                system(paste("rm", x))
            }
            if (gzipResults) {
                system(paste("gzip", filename.psl))
                filename.psl <- paste(filename.psl, "gz", sep = ".")
            }
            filename.psl
        }
    }
    return(unlist(filenames))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
