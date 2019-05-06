## ---- message=FALSE------------------------------------------------------
library(hiReadsProcessor)
## added to avoid package build failure on windows machines ##
if(.Platform$OS.type == "windows") { register(SerialParam()) }

## ------------------------------------------------------------------------
runData <- system.file("extdata/FLX_sample_run/", package = "hiReadsProcessor")
list.files(runData, recursive  = TRUE)

## ------------------------------------------------------------------------
seqProps <- read.SeqFolder(runData, seqfilePattern = ".+fna.gz$")
seqProps

## ------------------------------------------------------------------------
seqProps$sectors$"1"$samples

## ------------------------------------------------------------------------
seqProps <- findBarcodes(seqProps, sector = "all", showStats = TRUE)
seqProps
seqProps$sectors

## ------------------------------------------------------------------------
load(file.path(system.file("data", package = "hiReadsProcessor"),
               "FLX_seqProps.RData"))

## ---- eval=FALSE---------------------------------------------------------
#  seqProps <- findPrimers(seqProps, showStats=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  seqProps <- findLTRs(seqProps, showStats=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  seqProps <- findVector(seqProps, showStats = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  seqProps <- findLinkers(seqProps, showStats=TRUE, doRC=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  seqProps <- findIntegrations(seqProps,
#                               genomeIndices = c("hg18" = "/usr/local/genomeIndexes/hg18.noRandom.2bit"),
#                               numServers = 2)

## ------------------------------------------------------------------------
sampleSummary(seqProps)

## ------------------------------------------------------------------------
sessionInfo()

