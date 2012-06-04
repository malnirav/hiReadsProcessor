\name{clusterSites}
\alias{clusterSites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clusterSites(posID = NULL, value = NULL, grouping = NULL, psl.rd = NULL, weight = NULL, windowSize = 5, byQuartile = FALSE, quartile = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{posID}{
%%     ~~Describe \code{posID} here~~
}
  \item{value}{
%%     ~~Describe \code{value} here~~
}
  \item{grouping}{
%%     ~~Describe \code{grouping} here~~
}
  \item{psl.rd}{
%%     ~~Describe \code{psl.rd} here~~
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{windowSize}{
%%     ~~Describe \code{windowSize} here~~
}
  \item{byQuartile}{
%%     ~~Describe \code{byQuartile} here~~
}
  \item{quartile}{
%%     ~~Describe \code{quartile} here~~
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
function (posID = NULL, value = NULL, grouping = NULL, psl.rd = NULL, 
    weight = NULL, windowSize = 5, byQuartile = FALSE, quartile = 0.7) 
{
    if (is.null(psl.rd)) {
        stopifnot(!is.null(posID))
        stopifnot(!is.null(value))
    }
    else {
        if (!"Position" \%in\% colnames(psl.rd)) {
            stop("The object supplied in psl.rd parameter does not have Position column in it. Did you run getIntegrationSites() on it?")
        }
        posIDs <- paste(space(psl.rd), psl.rd$strand, sep = "")
        values <- psl.rd$Position
        if (is.null(weight)) {
            weight <- suppressWarnings(as.numeric(sub(".+counts=(.+)", 
                "\\1", psl.rd$qName)))
            if (all(is.na(weight))) {
                weight <- NULL
            }
            else {
                weight[is.na(weight)] <- 1
            }
        }
        grouping <- if (is.null(grouping)) {
            ""
        }
        else {
            grouping
        }
        clusters <- clusterSites(posIDs, values, grouping = grouping, 
            windowSize = windowSize, weight = weight, byQuartile = byQuartile, 
            quartile = quartile)
        message("Adding clustered data back to psl.rd.")
        clusteredValues <- with(clusters, split(clusteredValue, 
            paste(posID, value, grouping, sep = "")))
        psl.rd$clusteredPosition <- as.numeric(clusteredValues[paste(posIDs, 
            values, grouping, sep = "")])
        clusteredValueFreq <- with(clusters, split(clusteredValue.freq, 
            paste(posID, value, grouping, sep = "")))
        psl.rd$clonecount <- as.numeric(clusteredValueFreq[paste(posIDs, 
            values, grouping, sep = "")])
        rm("clusteredValueFreq", "clusteredValues", "clusters")
        cleanit <- gc()
        message("Picking best scoring hit to represent a cluster.")
        posIDs <- paste(space(psl.rd), psl.rd$strand, psl.rd$clusteredPosition, 
            sep = "")
        bestScore <- tapply(psl.rd$score, posIDs, max)
        isBest <- psl.rd$score == bestScore[posIDs]
        tocheck <- which(isBest)
        res <- tapply(tocheck, names(tocheck), "[[", 1)
        psl.rd$clusterTopHit <- FALSE
        psl.rd$clusterTopHit[res] <- TRUE
        message("Cleaning up!")
        rm("isBest", "bestScore", "posIDs", "tocheck", "res")
        cleanit <- gc()
        return(psl.rd)
    }
    groups <- if (is.null(grouping)) {
        ""
    }
    else {
        grouping
    }
    weight2 <- if (is.null(weight)) {
        1
    }
    else {
        weight
    }
    sites <- count(orderBy(~posID2 + value, data.frame(posID, 
        value, grouping = groups, weight = weight2, posID2 = paste(groups, 
            posID, sep = ""), stringsAsFactors = FALSE)), wt_var = "weight")
    rm("groups", "weight2")
    if (byQuartile) {
        message("Clustering by quartile: ", quartile)
        sites <- orderBy(~posID2 + value - freq, sites)
        quartiles <- with(sites, tapply(freq, posID2, quantile, 
            probs = quartile, names = FALSE))
        sites$belowQuartile <- with(sites, freq < quartiles[posID2])
        rm(quartiles)
        if (any(sites$belowQuartile)) {
            pos.be <- with(subset(sites, belowQuartile, drop = TRUE), 
                GRanges(IRanges(start = value, width = 1), seqnames = posID2, 
                  freq = freq))
            pos.ab <- with(subset(sites, !belowQuartile, drop = TRUE), 
                GRanges(IRanges(start = value, width = 1), seqnames = posID2, 
                  freq = freq))
            pos.overlap <- as.data.frame(as.matrix(findOverlaps(pos.be, 
                pos.ab, maxgap = windowSize, ignore.strand = TRUE)))
            pos.overlap$freq <- values(pos.ab[pos.overlap[, "subjectHits"]])$freq
            maxs <- with(pos.overlap, tapply(freq, as.character(queryHits), 
                max))
            pos.overlap$isMax <- with(pos.overlap, freq == maxs[as.character(queryHits)])
            rm(maxs)
            pos.overlap <- subset(pos.overlap, isMax, drop = TRUE)
            counts <- xtabs(isMax ~ queryHits, pos.overlap)
            if (length(table(counts)) > 1) {
                toFix <- as.numeric(names(which(counts > 1)))
                rows <- pos.overlap$queryHits \%in\% toFix
                pos.overlap$aboveQuartileValue <- pos.overlap$belowQuartileValue <- pos.overlap$valueDiff <- 0
                pos.overlap$aboveQuartileValue[rows] <- start(pos.ab[pos.overlap[rows, 
                  "subjectHits"]])
                pos.overlap$belowQuartileValue[rows] <- start(pos.be[pos.overlap[rows, 
                  "queryHits"]])
                pos.overlap$valueDiff[rows] <- with(pos.overlap[rows, 
                  ], abs(aboveQuartileValue - belowQuartileValue))
                mins <- with(pos.overlap[rows, ], tapply(valueDiff, 
                  as.character(queryHits), min))
                pos.overlap$isClosest <- TRUE
                pos.overlap$isClosest[rows] <- with(pos.overlap[rows, 
                  ], valueDiff == mins[as.character(queryHits)])
                pos.overlap <- subset(pos.overlap, isMax & isClosest, 
                  drop = TRUE)
                rm("counts", "mins")
            }
            pos.overlap$clusteredValue <- start(pos.ab[pos.overlap[, 
                "subjectHits"]])
            pos.overlap$posID2 <- as.character(seqnames(pos.be[pos.overlap[, 
                "queryHits"]]))
            rows <- which(!1:length(pos.be) \%in\% pos.overlap$query)
            loners <- pos.be[rows]
            if (length(loners) > 0) {
                times.rep <- values(loners)[["freq"]]
                res <- clusterSites(rep(as.character(seqnames(loners)), 
                  times = times.rep), rep(start(loners), times = times.rep), 
                  byQuartile = FALSE)
            }
            pos.overlap <- rbind(pos.overlap[, c("queryHits", 
                "clusteredValue")], data.frame(queryHits = rows, 
                clusteredValue = as.numeric(res$clusteredValue)))
            sites$clusteredValue <- sites$value
            sites$clusteredValue[sites$belowQuartile][pos.overlap[, 
                "queryHits"]] <- pos.overlap$clusteredValue
            stopifnot(any(!is.na(sites$clusteredValue)))
        }
        else {
            message("No sites found below defined quartile. Try to increase the quartile or use standard clustering, byQuartile=FALSE.")
        }
    }
    else {
        message("Clustering by minimum overlap.")
        sites.rl <- with(sites, RangedData(space = posID2, IRanges(start = value, 
            width = 1), freq))
        res <- as.data.frame(as.matrix(findOverlaps(sites.rl, 
            ignoreSelf = TRUE, ignoreRedundant = FALSE, select = "all", 
            maxgap = windowSize)))
        res$q.val <- start(sites.rl)[res$queryHits]
        res$s.val <- start(sites.rl)[res$subjectHits]
        res$q.freq <- sites.rl$freq[res$queryHits]
        res$s.freq <- sites.rl$freq[res$subjectHits]
        res$dist <- with(res, abs(q.val - s.val))
        stopifnot(!any(res$dist > 5))
        res$val <- with(res, ifelse(q.freq == s.freq, ifelse(q.val < 
            s.val, q.val, s.val), ifelse(q.freq >= s.freq, q.val, 
            s.val)))
        res$maxFreq <- with(res, pmax(q.freq, s.freq))
        maxes <- with(res, tapply(maxFreq, queryHits, max))
        res$ismaxFreq <- with(res, maxFreq == maxes[as.character(queryHits)])
        res <- orderBy(~-queryHits, res)
        hash.df <- unique(subset(res, ismaxFreq)[, c("queryHits", 
            "val")])
        clustered <- structure(as.numeric(hash.df$val), names = as.character(hash.df$queryHits))
        rm(hash.df)
        sites$clusteredValue <- sites$value
        sites$clusteredValue[as.numeric(names(clustered))] <- as.numeric(clustered)
        rm("clustered", "res")
        cleanit <- gc()
    }
    message("\t - Adding clustered value frequencies.")
    counts <- count(sites[, -grep("value", names(sites), fixed = TRUE)], 
        c("posID2", "clusteredValue"), wt_var = "freq")
    names(counts)[grep("freq", names(counts), fixed = TRUE)] <- "clusteredValue.freq"
    sites <- merge(sites, counts)
    if (byQuartile) {
        sites <- sites[, c("posID", "value", "freq", "clusteredValue", 
            "clusteredValue.freq", "grouping", "weight")]
    }
    sites$posID2 <- NULL
    if (is.null(grouping)) {
        sites$grouping <- NULL
    }
    if (is.null(weight)) {
        sites$weight <- NULL
    }
    return(sites)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
