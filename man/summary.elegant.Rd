\name{summary.elegant}
\alias{summary.elegant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.elegant(sampleInfo, samplenames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfo}{
%%     ~~Describe \code{sampleInfo} here~~
}
  \item{samplenames}{
%%     ~~Describe \code{samplenames} here~~
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
function (sampleInfo, samplenames = NULL) 
{
    stopifnot(class(sampleInfo) == "SimpleList")
    cat("Names of all sectors processed:", paste(names(sampleInfo$sectors), 
        collapse = ","), "\n")
    for (sector in names(sampleInfo$sectors)) {
        cat(rep("*", 50), "\n")
        cat("Sector:", sector, "\n")
        sampleList <- as.character(extractFeature(sampleInfo, 
            sector = sector, feature = "samplename")[[sector]])
        if (!is.null(samplenames)) {
            sampleList <- sampleList[sampleList \%in\% samplenames]
        }
        for (sample.i in sampleList) {
            cat("\t", sample.i, ":\n\t")
            metadatalist <- as.character(unlist(strsplit(extractFeature(sampleInfo, 
                sector = sector, samplename = sample.i, feature = "metadata")[[sector]], 
                ", ")))
            totalDecoded <- extractFeature(sampleInfo, sector = sector, 
                samplename = sample.i, feature = "decoded")[[sector]]
            if (length(totalDecoded) > 0) {
                totalDecoded <- length(totalDecoded[[1]])
                meta.processed <- grep("ed$", metadatalist, value = TRUE)
                meta.sites <- c("psl", "sites")
                for (metadata.i in metadatalist) {
                  cat("| ", metadata.i, ": ")
                  feature.i <- extractFeature(sampleInfo, sector = sector, 
                    samplename = sample.i, feature = metadata.i)[[sector]][[1]]
                  if (metadata.i \%in\% meta.processed) {
                    cat(length(feature.i), "(", round(100 * (length(feature.i)/totalDecoded), 
                      2), "\%) ", sep = "", "\n")
                  }
                  else if (metadata.i \%in\% meta.sites) {
                    cat(length(unique(feature.i$qName)), "(", 
                      round(100 * (length(unique(feature.i$qName))/totalDecoded), 
                        2), "\%) ", sep = "", "\n")
                  }
                  else {
                    cat(feature.i, " ")
                  }
                }
            }
            else {
                cat("\t None decoded. \n")
            }
            cat("\n", rep("~", 40), "\n")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
