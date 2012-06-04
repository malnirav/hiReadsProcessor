\name{otuSites}
\alias{otuSites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
otuSites(posID = NULL, readID = NULL, grouping = NULL, psl.rd = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{posID}{
%%     ~~Describe \code{posID} here~~
}
  \item{readID}{
%%     ~~Describe \code{readID} here~~
}
  \item{grouping}{
%%     ~~Describe \code{grouping} here~~
}
  \item{psl.rd}{
%%     ~~Describe \code{psl.rd} here~~
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
function (posID = NULL, readID = NULL, grouping = NULL, psl.rd = NULL) 
{
    if (is.null(psl.rd)) {
        stopifnot(!is.null(posID))
        stopifnot(!is.null(readID))
    }
    else {
        if (!"clusterTopHit" \%in\% colnames(psl.rd)) {
            stop("The object supplied in psl.rd parameter does not have Position column in it. Did you run clusterSites() on it?")
        }
        good.rows <- psl.rd$clusterTopHit
        posIDs <- paste(space(psl.rd), psl.rd$strand, psl.rd$clusteredPosition, 
            sep = "")
        readIDs <- psl.rd$qName
        grouping <- if (is.null(grouping)) {
            rep("A", nrow(psl.rd))
        }
        else {
            grouping
        }
        otus <- otuSites(posIDs[good.rows], readIDs[good.rows], 
            grouping[good.rows])
        message("Adding otuIDs back to psl.rd.")
        otuIDs <- with(otus, split(otuID, paste(posID, readID, 
            grouping)))
        psl.rd$otuIDs <- NA
        psl.rd$otuIDs[good.rows] <- as.numeric(otuIDs[paste(posIDs[good.rows], 
            readIDs[good.rows], grouping[good.rows])])
        rm("otus", "otuIDs", "posIDs", "readIDs", "grouping")
        cleanit <- gc()
        return(psl.rd)
    }
    groups <- if (is.null(grouping)) {
        "A"
    }
    else {
        grouping
    }
    sites <- data.frame(posID, readID, grouping = groups, stringsAsFactors = FALSE)
    rm(groups)
    reads <- ddply(sites, .(grouping, readID), summarise, posIDs = paste(sort(unique(posID)), 
        collapse = ","), counts = length(unique(posID)), .parallel = FALSE)
    reads <- orderBy(~grouping + posIDs, reads)
    reads$otuID <- unlist(lapply(lapply(with(reads, split(posIDs, 
        grouping)), as.factor), as.numeric))
    reads$newotuID <- reads$otuID
    singles <- reads$counts == 1
    toCheck <- with(reads[singles, ], split(posIDs, grouping))
    toCheck.ids <- with(reads[singles, ], split(otuID, grouping))
    toCheck <- sapply(names(toCheck), function(x) {
        names(toCheck[[x]]) <- toCheck.ids[[x]]
        toCheck[[x]]
    })
    rm(toCheck.ids)
    allposIDs <- with(reads[!singles, ], split(posIDs, grouping))
    for (f in intersect(names(toCheck), names(allposIDs))) {
        query <- structure(paste(toCheck[[f]], ",", sep = ""), 
            names = names(toCheck[[f]]))
        subject <- paste(allposIDs[[f]], ",", sep = "")
        res <- sapply(query, grep, x = subject, fixed = TRUE)
        res <- res[sapply(res, length) > 0]
        res <- structure(unlist(res, use.names = F), names = rep(names(res), 
            sapply(res, length)))
        reads[!singles & reads$grouping == f, "newotuID"][as.numeric(res)] <- as.numeric(names(res))
    }
    ots.ids <- with(reads, split(newotuID, readID))
    sites$otuID <- as.numeric(unlist(ots.ids[sites$readID]))
    stopifnot(any(!is.na(sites$otuID)))
    rm(reads)
    cleanit <- gc()
    if (is.null(grouping)) {
        sites$grouping <- NULL
    }
    return(sites)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
