#' Functions to process LM-PCR reads from 454/Solexa data.
#'
#' hiReadsProcessor contains set of functions which allow users to process LM-PCR sequence data coming out of the 454/Solexa sequencer. Given an excel file containing parameters for demultiplexing and sample metadata, the functions automate trimming of adaptors and identification of the genomic product. In addition, if IntSites MySQL database is setup, the sequence attrition is loaded into respective tables for post processing setup and analysis.
#'
#' @import Biostrings GenomicRanges foreach iterators RMySQL xlsx plyr Rsubread ShortRead
#' @docType package
#' @name hiReadsProcessor
NULL

#' Read contents of a sequencing folder and make a SimpleList object
#'
#' Given a sequencing folder path, sample information file path, and sequence file extension pattern, the function returns a list of variables required to process the data. The function also calls \code{\link{read.sampleInfo}} which reads in sample processing metadata and formats it if needed.
#'
#' @param sequencingFolderPath full or relative path to the sequencing folder
#' @param sampleInfoFilePath full or relative path to the sample information file, which holds samples to quadrant/lane associations along with other metadata required to trim sequences or process it. Default to NULL, where the function tries to find xls or tab deliminated txt file in the sequencing folder which sounds similar to 'sampleinfo' and present you with choices of file to select from.
#' @param seqfilePattern regex to describe sequence file endings. See examples. Default is "\\.TCA.454Reads.fna$".
#' @param interactive whether to prompt each time the function encounters an issue or use the defaults. Default is TRUE.
#'
#' @return a SimpleList list which is used by other functions to process and decode the data.
#'
#' @note 
#' \itemize{
#'   \item One must make sure that each sequencing file has sector name/number prefixed at the beginning, else \code{\link{decodeByBarcode}} will fail trying to find the filename.
#'   \item For paired end Illumina runs, make sure the filenames include R1, R2, and I1 somewhere in the name denoting pair1, pair2, and index/barcode reads, respectively. 
#' }
#'
#' @seealso \code{\link{read.sampleInfo}}, \code{\link{decodeByBarcode}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples  
#' #read.SeqFolder("~/Downloads/454Runs/2011_07_22", seqfilePattern=".+fna$")
#' #read.SeqFolder(".", seqfilePattern=".+fastq$")
#'
read.SeqFolder <- function(sequencingFolderPath=NULL, sampleInfoFilePath=NULL, 
                           seqfilePattern="\\.TCA.454Reads.fna$", interactive=TRUE) {
  if(is.null(sequencingFolderPath)) {
    stop("No Sequencing Folder Path provided.")
  }
  
  ## get the sequencingFolderPath right!
  sequencingFolderPath <- normalizePath(sequencingFolderPath, mustWork=TRUE)
  
  seqFilePaths <- list.files(path=sequencingFolderPath, recursive=TRUE, 
                             full.names=TRUE, pattern=seqfilePattern)
  if(length(seqFilePaths)==0) {
    stop("No files found in the folder ", sequencingFolderPath,
         "matching following pattern: ", seqfilePattern)
  }
  
  ## read the sample info file
  if(is.null(sampleInfoFilePath)) {
    possibleFiles <- list.files(path=sequencingFolderPath, 
                                recursive=TRUE, full.names=TRUE, 
                                pattern=".*sampleinfo.+", ignore.case=TRUE)
    if (length(possibleFiles)==0) {
      stop("No sample information file found in folder: ", sequencingFolderPath)
    } else {
      if(interactive & length(possibleFiles)>1) {
        message("Please choose a sample information file to read the meta data from:\n",
                paste(1:length(possibleFiles), possibleFiles, sep=": ", collapse="\n"))
        choice <- scan(what=integer(0), n=1, quiet=TRUE, multi.line=FALSE)
      } else {
        choice <- 1            
      }
      message("Choosing ", possibleFiles[choice], " as sample information file.")
    }
    sampleInfoFilePath <- possibleFiles[choice]
  }
  sampleInfo <- read.sampleInfo(sampleInfoFilePath, interactive=interactive)
  
  ## do a quick test of filenames if any samples are from paired end illumina
  if(any(sampleInfo$pairedend)) {    
    sectors <- unique(sampleInfo$sector[sampleInfo$pairedend])
    for(sector in sectors) {
      vals <- grep(gsub("I1|R1|R2","",sector), seqFilePaths, value=TRUE)
      if(length(vals)!=3) {
        stop("Sector ",sector," is missing one of the files: R1, R2, or I1.")
      }
    }
  }
  
  if(length(sampleInfo)!=length(unique(gsub("I1|R1|R2","",seqFilePaths)))) {
    warning("Number of sectors (", length(sampleInfo),
            ") in sample information file does not match number of sector files (", 
            length(unique(gsub(seqfilePattern,'',seqFilePaths))), 
            ") found in the folder.")
  }
  
  return(SimpleList("sequencingFolderPath"=sequencingFolderPath, 
                    "seqFilePaths"=seqFilePaths, "seqfilePattern"=seqfilePattern, 
                    "sampleInfoFilePath"=sampleInfoFilePath, 
                    "sectors"=sampleInfo, "callHistory"=match.call()))
}

#' Read a sample information file and format appropriate metadata.
#'
#' Given a sample information file, the function checks if it includes required information to process samples present on each sector/quadrant/region/lane. The function also adds other columns required for processing with default values if not already defined ahead of time.
#'
#' @param sampleInfoFilePath full or relative path to the sample information file, which holds samples to quadrant/lane associations along with other metadata required to trim sequences or process it. 
#' @param splitBySector split the data frame into a list by sector column. Default is TRUE.
#' @param interactive whether to prompt each time the function encounters an issue, or use the defaults. Default is TRUE.
#'
#' @details
#' \itemize{
#'  \item Required Column Description:
#'    \itemize{
#'  	\item sector => region/quadrant/lane of the sequencing plate the sample comes from. If files have been split by samples apriori, then the filename associated per sample without the extension. If this is a filename, then be sure to enable 'alreadyDecoded' parameter in \code{\link{findBarcodes}} or \code{\link{decodeByBarcode}}, since contents of this column is pasted together with 'seqfilePattern' parameter in \code{\link{read.SeqFolder}} to find the appropriate file needed. For paired end data, this is basename of the FASTA/Q file holding the sample data from the LTR side. For example, files such as Lib3_L001_R2_001.fastq.gz or Lib3_L001_R2_001.fastq would be Lib3_L001_R2_001, and consequently Lib3_L001_R1_001 would be used as the second pair!
#'  	\item barcode => unique 4-12bp DNA sequence which identifies the sample. If providing filename as sector, then leave this blank since it is assumed that the data is already demultiplexed.
#'  	\item primerltrsequence => DNA sequence of the viral LTR primer with/without the viral LTR sequence following the primer landing site. If already trimmed, then mark this as SKIP.
#'  	\item samplename => Name of the sample associated with the barcode
#'  	\item sampledescription => Detailed description of the sample
#'  	\item gender => sex of the sample: male or female or NA
#'  	\item species => species of the sample: homo sapien, mus musculus, etc.
#'  	\item freeze => UCSC freeze to which the sample should be aligned to.
#'  	\item linkersequence => DNA sequence of the linker adaptor following the genomic sequence. If already trimmed, then mark this as SKIP.
#'  	\item restrictionenzyme => Restriction enzyme used for digestion and sample recovery. Can also be one of: Fragmentase or Sonication!
#' 		}
#'  \item Metadata Parameter Column Description:
#'   \itemize{
#'  	\item ltrbitsequence => DNA sequence of the viral LTR following the primer landing site. Default is last 7bps of the primerltrsequence.
#'  	\item ltrbitidentity => percent of LTR bit sequence to match during the alignment. Default is 1.
#'  	\item primerltridentity => percent of primer to match during the alignment. Default is .85
#'  	\item linkeridentity => percent of linker sequence to match during the alignment. Default is 0.55. Only applies to non-primerID/random tag based linker search.
#'  	\item primeridinlinker => whether the linker adaptor used has primerID/random tag in it? Default is FALSE.
#'  	\item primeridinlinkeridentity1 => percent of sequence to match before the random tag. Default is 0.75. Only applies to primerID/random tag based linker search and when primeridinlinker is TRUE.
#'  	\item primeridinlinkeridentity2 => percent of sequence to match after the random tag. Default is 0.50. Only applies to primerID/random tag based linker search and when primeridinlinker is TRUE.
#'  	\item celltype => celltype information associated with the sample
#'  	\item user => name of the user who prepared or processed the sample
#'  	\item pairedend => is the data paired end? Default is FALSE.
#' 		}
#'  \item Processing Parameter Column Description:
#'   \itemize{
#'  	\item startwithin => upper bound limit of where the alignment should start within the query. Default is 3.
#'  	\item alignratiothreshold => cuttoff for (alignment span/read length). Default is 0.7.
#'  	\item genomicpercentidentity => cuttoff for (1-(misMatches/matches)). Default is 0.98.
#'  	\item clustersiteswithin => cluster integration sites within a defined window size based on frequency which corrects for any sequencing errors. Default is 5.
#'  	\item keepmultihits => whether to keep sequences/reads that return multiple best hits, aka ambiguous locations. 
#'  	\item processingdate => the date of processing 
#' 		}
#' }
#'
#' @return if splitBySector=TRUE, then an object of SimpleList named by quadrant/lane information defined in sampleInfo file, else a dataframe.
#'
#' @seealso \code{\link{read.SeqFolder}}, \code{\link{decodeByBarcode}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples  
#' #read.sampleInfo("~/sampleInfo.xls")
#'
read.sampleInfo <- function(sampleInfoPath=NULL, splitBySector=TRUE, 
                            interactive=TRUE) {
  ## read file and make sampleInfo object with sample to metadata associations ##
  if(is.null(sampleInfoPath)) {
    stop("No sample information file path provided.")
  }
  
  sampleInfoPath <- normalizePath(sampleInfoPath, mustWork=TRUE)
  
  requiredCols <- c('sector', 'barcode', 'primerltrsequence', 'samplename', 
                    'sampledescription', 'gender', 'species', 'freeze', 
                    'linkersequence', 'restrictionenzyme')
  
  metaDataCols <- c('ltrbitsequence'='', 'ltrbitidentity'=1, 
                    'primerltridentity'=.85, 'linkeridentity'=.55, 
                    'primeridinlinker'=FALSE, 'primeridinlinkeridentity1'=.75, 
                    'primeridinlinkeridentity2'=.50, 'celltype'='', 
                    'user'=Sys.getenv("USER"), 'startwithin'=3, 
                    'alignratiothreshold'=.7, 'clustersiteswithin'=5, 
                    'keepmultihits'=TRUE , 'genomicpercentidentity'=.98, 
                    'processingdate'=format(Sys.time(), "%Y-%m-%d "), 
                    'pairedend'=FALSE)
  
  if(grepl('.xls.?$', sampleInfoPath)) {
    sampleInfo <- unique(read.xlsx(sampleInfoPath, 
                                   sheetIndex=1, stringsAsFactors=FALSE))
  } else {
    sampleInfo <- unique(read.delim(sampleInfoPath, stringsAsFactors=FALSE))
  }
  names(sampleInfo) <- tolower(gsub("\\.|-|_", "", names(sampleInfo)))
  
  # check for required columns
  ColsNotThere <- !requiredCols %in% names(sampleInfo)
  if (any(ColsNotThere)) {
    absentCols <- requiredCols[ColsNotThere]
    stop("Following required column(s) is absent from the Sample Info file: ",
         paste(absentCols,sep="", collapse=", "))
  }
  
  # add missing meta data columns
  metaColsNotThere <- !names(metaDataCols) %in% names(sampleInfo)
  if(any(metaColsNotThere)) {
    sampleInfo <- cbind(sampleInfo,
                        as.data.frame(t(metaDataCols[metaColsNotThere]),
                                      stringsAsFactors = FALSE))        
  }
  
  # do some formatting to avoid later hassels!
  for(column in c('sector', 'barcode', 'primerltrsequence', 'ltrbitsequence', 
                  'samplename', 'linkersequence', 'restrictionenzyme')) {
    sampleInfo[,column] <- gsub(" ", "", sampleInfo[,column])
    if(column %in% c('barcode', 'primerltrsequence', 'ltrbitsequence', 
                     'linkersequence', 'restrictionenzyme')) {
      sampleInfo[,column] <- toupper(sampleInfo[,column])
    }
  }
  
  for(column in c('pairedend', 'keepmultihits', 'primeridinlinker')) {
    sampleInfo[,column] <- as.logical(sampleInfo[,column])
  }
  
  # confirm ltr bit is correct
  ltrbitTest <- sampleInfo$primerltrsequence=="SKIP"
  if(any(ltrbitTest)) { ## add SKIP to ltrbit as well if primerltrsequence has been trimmed
    tofix <- which(ltrbitTest)
    message("adding SKIP to ltrbitsequence to ",length(tofix),
            " sample since primerltrsequence has been trimmed.")
    sampleInfo$ltrbitsequence[tofix] <- "SKIP"
  }
  
  ltrbitTest <- nchar(sampleInfo$ltrbitsequence)==0 | sampleInfo$ltrbitsequence==""
  if(any(ltrbitTest)) {
    tofix <- which(ltrbitTest)
    if(interactive) {
      message("LTR bit not found for ",length(tofix),
              " samples. Use last 7 bases of the LTR primer as the LTR bit? (y/n)")
      choice <- scan(what=character(0), n=1, quiet=TRUE, multi.line=FALSE)
    } else {
      message("LTR bit not found for ",length(tofix),
              " samples. Using last 7 bases of the LTR primer as the LTR bit.")
      choice <- "y"
    }
    
    if(tolower(choice)=="y") {
      sampleInfo$ltrbitsequence <- substr(sampleInfo$primerltrsequence,
                                          nchar(sampleInfo$primerltrsequence)-6, 
                                          nchar(sampleInfo$primerltrsequence))
      sampleInfo$primerltrsequence <- substr(sampleInfo$primerltrsequence, 1,
                                             nchar(sampleInfo$primerltrsequence)-7)
    } else {
      warning("No LTR bit sequence found for following samples: ",
              paste(sampleInfo$samplename[tofix], sep="", collapse=", "),
              immediate.=TRUE)
    }        
  }
  
  # check if samplenames are up to the expectations
  samplenametest <- nchar(sampleInfo$samplename)==0 | sampleInfo$samplename==""
  if(any(samplenametest)) {
    stop("No sample names found in following rows of the sample information file ",
         sampleInfoPath, " : ", paste(which(samplenametest), sep="", collapse=", "))
  }
  
  # check for sectors and their usage
  sectortest <- nchar(sampleInfo$sector)==0 | sampleInfo$sector=="" | 
                is.na(sampleInfo$sector)
  if(any(sectortest)) {
    tofix <- which(sectortest)
    if(interactive) {
      message("Sector information not found for ", length(tofix),
              " samples. Which sector are they from? (1,2,4,etc)")
      choice <- scan(what=character(0), quiet=TRUE, multi.line=FALSE)
    } else {
      message("Sector information not found for ", length(tofix),
              " samples. Assuming they are from sector 1.")
      choice <- "1"
    }
    if(length(choice)>0) {
      sampleInfo$sector[tofix] <- unlist(strsplit(choice,","))
    } else {
      stop("No Sector information found for following samples: ",
           paste(sampleInfo$samplename[tofix], sep="", collapse=", "))
    }
  }
  ## excel sometimes converts integers to doubles...make sure to remove the trailing 0
  sampleInfo$sector <- gsub("\\.0$", "", as.character(sampleInfo$sector))
  
  sampleSectorTest <- table(paste(sampleInfo$samplename, sampleInfo$sector))
  if(any(sampleSectorTest>1)) {
    stop("Duplicate sample names found on same quadrant in the ",
         "sample information file ", sampleInfoPath, " : ", 
         paste(sampleSectorTest[sampleSectorTest>1], sep="",collapse=", "))
  }
  
  # prepare the sample info object!
  if(splitBySector) {            
    sampleInfo <- SimpleList(split(sampleInfo, sampleInfo$sector))
    for(sector in 1:length(sampleInfo)) { 
      sampleData <- SimpleList(split(sampleInfo[[sector]], 
                                     sampleInfo[[sector]]$samplename))
      for(sample.i in 1:length(sampleData)) { 
        sampleData[[sample.i]] <- SimpleList(as.list(sampleData[[sample.i]])) 
      }
      sampleInfo[[sector]] <- SimpleList("samples"=sampleData) 
    }
  }
  return(sampleInfo)
}

#' Removes duplicate sequences from DNAStringSet object.
#'
#' Given a DNAStringSet object, the function dereplicates reads and adds counts=X to the definition line to indicate replication. 
#'
#' @param dnaSet DNAStringSet object to dereplicate. 
#'
#' @return DNAStringSet object with names describing frequency of repeat.
#'
#' @seealso \code{\link{replicateReads}}, \code{\link{removeReadsWithNs}}, \code{\link{decodeByBarcode}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples  
#' #dereplicateReads(dnaSet)
#'
dereplicateReads <- function(dnaSet) {
  stopifnot(is(dnaSet,"DNAStringSet"))
  if(is.null(names(dnaSet))) {
    message("No names attribute found in dnaSet object...", 
            "using artifically generated names")
    names(dnaSet) <- paste("read", 1:length(dnaSet), sep="-")
  }
  dnaSet <- sort(dnaSet)
  ranks <- rank(dnaSet)
  counts <- table(ranks)
  isDuplicate <- duplicated(ranks)
  seqToRank <- data.frame(ranks, counts=as.numeric(counts[as.character(ranks)]),
                          row.names=names(dnaSet), stringsAsFactors = FALSE)
  seqToRank <- seqToRank[!isDuplicate,]
  dnaSet <- dnaSet[!isDuplicate]    
  names(dnaSet) <- paste0(names(dnaSet), 
                          "counts=", seqToRank[names(dnaSet),"counts"])
  return(dnaSet)
}

#' Replicate sequences from DNAStringSet object using counts identifier or vector
#'
#' Given a DNAStringSet object, the function replicates reads using counts=X marker at the end of definition line. 
#'
#' @param dnaSet DNAStringSet object to replicate. 
#' @param counts an integer or a numeric vector of length length(dnaSet) indicating how many times to repeat each sequence. Default is NULL, in which it uses counts=X notation from the definition line to replicate reads.
#'
#' @return DNAStringSet object.
#'
#' @seealso \code{\link{dereplicateReads}}, \code{\link{removeReadsWithNs}}, \code{\link{decodeByBarcode}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples  
#' #replicateReads(decoded)
#'
replicateReads <- function(dnaSet, counts=NULL) {
  stopifnot(is(dnaSet,"DNAStringSet"))
  if(is.null(counts)) {
    if(is.null(names(dnaSet))) {
      stop("No names attribute found in dnaSet object")
    }
    counts <- as.numeric(sub(".+counts=(\\d+)","\\1", names(dnaSet)))
    if(all(is.na(counts))) {
      stop("No counts=X marker found at the end of definition line or ",
           "names attribute in dnaSet object")
    }
  }
  if (length(counts)==1) {
    counts <- rep(counts, length(dnaSet))
  }
  
  ids <- unlist(sapply(counts, seq_len))
  deflines <- paste0(rep(sub("(.+)counts=.+", "\\1", 
                             names(dnaSet)), times=counts), ids)
  dnaSet <- rep(dnaSet, times=counts)
  names(dnaSet) <- deflines
  return(dnaSet)
}

#' Remove sequences with ambiguous nucleotides.
#'
#' Given a DNAStringSet object, the function removes any reads that has either repeating or total Ns which is greater than to maxNs threshold
#'
#' @param dnaSet DNAStringSet object to evaluate. 
#' @param maxNs integer value denoting the threshold of maximum allowed Ns. Default is 5.
#' @param consecutive boolean flag denoting whether Ns to filter is consecutive or total . Default is TRUE.
#'
#' @return DNAStringSet object.
#'
#' @seealso \code{\link{dereplicateReads}}, \code{\link{replicateReads}}, \code{\link{decodeByBarcode}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples  
#' #removeReadsWithNs(dnaSet)
#'
removeReadsWithNs <- function(dnaSet, maxNs=5, consecutive=TRUE) {
  if(consecutive) {
    good.row <- grepl(paste(rep("N",maxNs), collapse=""), dnaSet, fixed=TRUE)
  } else {
    res <- alphabetFrequency(dnaSet)
    good.row <- res[,"N"] <= maxNs
  }
  return(dnaSet[good.row])
}

#' Breaks an object into chunks of N size.
#'
#' Given a linear object, the function breaks query into chunks of N size where each chunk has a respective subject object filtered by seqnames/space present in the query chunk. This is a helper function used by functions in 'See Also' section where each chunk is sent to a parallel node for processing.
#'
#' @param x a linear object.
#' @param chunkSize number of rows to use per chunk of query. Default to length(sites.rd)/detectCores() or length(query)/getDoParWorkers() depending on parallel backend registered. 
#'
#' @return a list of object split into chunks.
#'
#' @seealso \code{\link{primerIDAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{pairwiseAlignSeqs}}
#'
#' @export
#'
#' @examples
#' x <- c("GAGGCTGTCACCGACAAGGTTCTGA", "AATAGCGTGGTGACAGCCCACATGC", 
#'        "GGTCTTCTAGGGAACCTACGCCACA", "TTTCCGGCGGCAGTCAGAGCCAAAG", 
#'        "TCCTGTCAACTCGTAGATCCAATCA", "ATGGTCACCTACACACAACGGCAAT", 
#'        "GTCAGGACACCTAATCACAAACGGA", "AGACGCAGGTTCAGGCTGGACTAGA", 
#'        "ATCGTTTCCGGAATTCGTGCACTGG", "CAATGCGGGCACACGCTCTTACAGT")
#' chunkize(DNAStringSet(x), 5)
chunkize <- function(x, chunkSize = NULL) {
  chunks <- breakInChunks(length(x), 
                          ifelse(!is.null(chunkSize),
                                 length(x)/chunkSize, 
                                 ifelse(!is.null(is.null(getDoParWorkers())), 
                                        length(x)/getDoParWorkers(), 
                                        length(x)/detectCores())))
  mapply(function(z, y) x[z:y], start(chunks), end(chunks), 
         SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

#' Split DNAStringSet object using first X number of bases defined by a vector.
#'
#' Given a character vector of barcodes/MID to sample association and a DNAStringSet object, the function splits/demultiplexes the DNAStringSet object by first few bases dictated by length of barcodes/MID supplied. This is an accessory function used by \code{\link{decodeByBarcode}} or \code{\link{findBarcodes}}
#'
#' @param barcodesSample a character vector of barcodes to sample name associations. Ex: c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2",...)
#' @param dnaSet DNAStringSet object to evaluate. 
#' @param trimFrom integer value serving as start point to trim the sequences from. This is calculated internally length barcode+1. Default is NULL.
#' @param showStats boolean flag denoting whether to show decoding statistics per sample & barcode. Default is FALSE.
#' @param returnUnmatched boolean flag denoting whether to return unmatched reads. Default is FALSE.
#' @param dereplicate return dereplicated sequences. Calls \code{\link{dereplicateReads}}, which appends counts=X to sequence names/deflines. Default is FALSE.
#'
#' @return DNAStringSet object split by sample name found in barcodesSample.
#'
#' @seealso \code{\link{decodeByBarcode}}, \code{\link{dereplicateReads}}, \code{\link{replicateReads}}
#'
#' @export
#'
#' @examples 
#' #splitByBarcode(c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2"), dnaSet, showStats=TRUE, returnUnmatched=TRUE)
#'
splitByBarcode <- function(barcodesSample, dnaSet, trimFrom=NULL, showStats=FALSE, 
                           returnUnmatched=FALSE, dereplicate=FALSE) {
  if(is.null(barcodesSample) | length(barcodesSample)==0) {
    stop("No barcodes to samples association vector provided in parameter ",
         "barcodesSample.")
  }
  
  if(is.null(dnaSet) | length(dnaSet)==0) {
    stop("No sequences provided in parameter dnaSet.")
  }
  
  if(is.null(names(dnaSet))) {
    stop("No names attribute found in dnaSet object")
  }
  
  message("Using following schema for barcode to sample associations")
  print(as.data.frame(barcodesSample))
  
  ## subset barcode string from the rest of sequence ##
  barcodelen <- unique(nchar(names(barcodesSample)))
  seqbarcodes <- substr(dnaSet,1,barcodelen)
  
  ## index rows that match your list of barcodes ##
  good.row <- seqbarcodes %in% names(barcodesSample)
  if(!any(good.row)) {
    stop("No matching barcoded sequences found on this quadrant.")
  }
  
  sampleNames <- barcodesSample[seqbarcodes[good.row]]    
  deflines <- sub("^(\\S+) .+$", "\\1", names(dnaSet)[good.row], perl=TRUE)
  
  ## if primer bases were utilized for tiebreaker, use the original length instead of modified for trimming.
  if(is.null(trimFrom)) {        
    trimFrom <- barcodelen+1
  }
  
  ## remove sequences with unknown barcode and trim barcode itself ##
  unmatched <- DNAStringSet(dnaSet[!good.row])
  dnaSet <- DNAStringSet(dnaSet[good.row], start=trimFrom)
  names(dnaSet) <- deflines
  
  if(showStats) {
    message("Number of Sequences with no matching barcode: ",
            as.numeric(table(good.row)['FALSE']))
    message("Number of Sequences decoded:")
    print(as.data.frame(table(sampleNames)))
  }
  
  dnaSet <- as.list(split(dnaSet, as.character(sampleNames)))
  
  if(dereplicate) {
    message("Dereplicating reads.")
    dnaSet <- lapply(dnaSet, dereplicateReads)
  }
  
  if(returnUnmatched) {
    dnaSet <- c(dnaSet, "unDecodedSeqs"=unmatched)
  }
  
  return(dnaSet)
}

#' Demultiplex reads by their barcodes
#'
#' Given a sample information object, the function reads in the raw fasta/fastq file, demultiplexes reads by their barcodes, and appends it back to the sampleInfo object. Calls \code{\link{splitByBarcode}} to perform the actual splitting of file by barcode sequences. If supplied with a character vector and reads themselves, the function behaves a bit differently. See the examples.
#'
#' @param sampleInfo sample information SimpleList object created using \code{\link{read.sampleInfo}}, which holds barcodes and sample names per sector/quadrant/lane or a character vector of barcodes to sample name associations. Ex: c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2",...)
#' @param sector If sampleInfo is a SimpleList object, then a numeric/character value or vector representing sector(s) in sampleInfo. Optionally if on high memory machine sector='all' will decode/demultiplex sequences from all sectors/quadrants. This option is ignored if sampleInfo is a character vector. Default is NULL. 
#' @param dnaSet DNAStringSet object containing sequences to be decoded or demultiplexed. Default is NULL. If sampleInfo is a SimpleList object, then reads are automatically extracted using \code{\link{read.seqsFromSector}} and parameters defined in sampleInfo object.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param returnUnmatched return unmatched sequences. Returns results as a list where x[["unDecodedSeqs"]] has culprits. Default is FALSE.
#' @param dereplicate return dereplicated sequences. Calls \code{\link{dereplicateReads}}, which appends counts=X to sequence names/deflines. Default is FALSE.
#' @param alreadyDecoded if reads have be already decoded and split into respective files per sample and 'seqfilePattern' parameter in \code{\link{read.SeqFolder}} is set to reading sample files and not the sector files, then set this to TRUE. Default is FALSE. Enabling this parameter skips the barcode detection step and loads the sequence file as is into the sampleInfo object. 
#'
#' @return If sampleInfo is an object of SimpleList then decoded sequences are appeneded to respective sample slots, else a named list of DNAStringSet object. If returnUnmatched=TRUE, then x[["unDecodedSeqs"]] has the unmatched sequences.
#'
#' @seealso \code{\link{splitByBarcode}}, \code{\link{dereplicateReads}}, \code{\link{replicateReads}}
#'
#' @export
#'
#' @aliases findBarcodes
#' 
#' @examples 
#' #findBarcodes(sampleInfo,showStats=TRUE)
#' #decodeByBarcode(sampleInfo,showStats=TRUE)
#' #decodeByBarcode(sampleInfo=c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2"), dnaSet, showStats=TRUE, returnUnmatched=TRUE)
#'
decodeByBarcode <- function(sampleInfo, sector=NULL, dnaSet=NULL, showStats=FALSE, 
                            returnUnmatched=FALSE, dereplicate=FALSE, 
                            alreadyDecoded=FALSE) {
  
  ## tried PDict()...and its slower than this awesome code! ##    
  if(is(sampleInfo,"SimpleList")) {
    if(is.null(sector)) {
      stop("No sector provided in parameter sector.")
    }
    
    sectors <- sector <- as.character(sector)
    if(length(sectors)==1 & tolower(sectors)=="all") {
      sectors <- names(sampleInfo$sectors)
    }
    
    if(any(!sectors %in% names(sampleInfo$sectors))) {
      stop("Following sectors not found in names(sampleInfo$sectors): ",
           sectors[!sectors %in% names(sampleInfo$sectors)])
    }
    
    for(sector in sectors) {
      ## check everything is cool with the provided barcodes first before reading the sequences! ##
      message("Decoding sector: ",sector) 
      isPaired <- any(as.logical(extractFeature(sampleInfo, sector, 
                                                feature='pairedend')[[sector]]))
      
      ## prepare a vector of barcode to sample associations ##
      sampleBarcodes <- toupper(extractFeature(sampleInfo, sector=sector,
                                               feature="barcode")[[sector]])
      barcodesSample <- structure(names(sampleBarcodes), 
                                  names=as.character(sampleBarcodes))
      
      if (length(table(nchar(as.character(sampleBarcodes))))>1) {
        stop("Multiple barcode lengths found.")
      }
      
      ## length of barcodes before any modifications done later if any! ##
      realbarcodelen <- unique(nchar(as.character(sampleBarcodes)))
      
      if (any(table(as.character(sampleBarcodes))>1)) {
        message("Duplicate barcode found on this sector.\n",
                "Please choose from one of the options below:\n",
                "\t1: Pick first few bases of primer for tiebreaker? ",
                "(This could be dangerous if the run has too many errors!)\n",
                "\t2: Use the last sample associated with the duplicate ",
                "as the primary sample?\n", 
                "\t3: Do not do anything.")
        
        choice <- scan(what=integer(0), n=1, quiet=TRUE, multi.line=FALSE)
        
        if(choice==1) {  
          message("Enter # of bases to use from primer:")
          howmany <- scan(what=integer(0), n=1, quiet=TRUE, multi.line=FALSE)
          samplePrimers <- extractFeature(sampleInfo,
                                          sector=sector,
                                          feature="primerltrsequence")[[sector]]
          samplePrimers <- toupper(samplePrimers)
          newBarcodes <- toupper(paste0(sampleBarcodes,
                                        substr(samplePrimers,1,howmany))) 
          counts <- table(newBarcodes)
          ## only take 1 to 1 associations!
          rows <- newBarcodes %in% names(which(counts==1)) 
          
          if (any(counts>1)) {            
            message("Tie breaking failed...",
                    "try choosing high number of bases from primer possibly? ",
                    "Here are the failed barcodes: ",
                    paste(names(which(counts>1)),collapse=", "))
            
            message("Ignore samples associated with those barcodes and ",
                    "continue processing? (y/n)")
            whatsup <- scan(what=character(0), n=1, quiet=TRUE, multi.line=FALSE)            
            if(whatsup=='n') {
              stop("Aborting processing due to ambiguous barcode ",
                   "association for samples: ",
                   paste(names(sampleBarcodes[!rows]),collapse=", "))
            } else {              
              message("Ignoring following samples due to duplicate barcodes: ",
                      paste(names(sampleBarcodes[!rows]),collapse=", "))
            }
          }
          
          barcodesSample <- structure(names(sampleBarcodes[rows]), 
                                      names=newBarcodes[rows])
          
        } else if(choice==2) {
          message("Overwriting duplicate samples associated with the same barcode...")
        } else {
          stop("Aborting due to duplicate barcode found on this sector")
        }
      }
      
      dnaSet <- read.seqsFromSector(sampleInfo, sector, isPaired)
      
      if(alreadyDecoded) {
        if(length(barcodesSample)>1) {
          stop("alreadyDecoded parameter is set to TRUE. There shouldn't be more ",
               "than one sample associated to a sequence file.")
        }
        
        if(is.list(dnaSet)) {
          dnaSet <- sapply(dnaSet, function(x) {
            names(x) <- sub("^\\S+-(\\S+) .+$", "\\1", names(x), perl=TRUE)  
            x
          })          
        } else {
          names(dnaSet) <- sub("^\\S+-(\\S+) .+$", "\\1", names(dnaSet), perl=TRUE) 
        }        
                
        if(isPaired) {
          ## no need to store barcode/index reads if alreadyDecoded!
          dnaSet <- dnaSet[c("pair1","pair2")]
          dnaSet <- sapply(dnaSet, function(x) {
            as.list(split(x, rep(as.character(barcodesSample), length(x))))
          })
          dnaSet <- mapply(function(x,y) list("pair1"=x, "pair2"=y), 
                           dnaSet[["pair1"]], dnaSet[["pair2"]], SIMPLIFY=FALSE)
        } else {
          dnaSet <- as.list(split(dnaSet,
                                  rep(as.character(barcodesSample), length(dnaSet))))
        }
      } else {
        if(isPaired) {
          bc <- splitByBarcode(barcodesSample, dnaSet[["barcode"]], 
                               trimFrom=realbarcodelen+1, showStats=showStats, 
                               returnUnmatched=returnUnmatched, 
                               dereplicate=dereplicate)
          p1 <- sapply(bc, function(x) dnaSet[['pair1']][names(x)])
          p2 <- sapply(bc, function(x) dnaSet[['pair2']][names(x)])
          stopifnot(identical(sapply(bc,length), sapply(p1,length)))
          stopifnot(identical(sapply(bc,length), sapply(p2,length)))
          dnaSet <- mapply(function(x,y) list("pair1"=x, "pair2"=y), p1, p2, 
                           SIMPLIFY=FALSE)
          rm("bc","p1","p2")
        } else {
          dnaSet <- splitByBarcode(barcodesSample, dnaSet, 
                                   trimFrom=realbarcodelen+1, showStats=showStats, 
                                   returnUnmatched=returnUnmatched, 
                                   dereplicate=dereplicate)
        }
      }
      
      for(samplename in names(dnaSet)) {
        if(samplename=="unDecodedSeqs") {                                        
          metadata(sampleInfo$sectors[[sector]]) <- 
            append(metadata(sampleInfo$sectors[[sector]]), 
                   list("unDecodedSeqs"=dnaSet[[samplename]]))            
        } else {
          sampleInfo$sectors[[sector]]$samples[[samplename]]$decoded <- 
            dnaSet[[samplename]]
        }
      }
      metadata(sampleInfo$sectors[[sector]]) <- 
        append(metadata(sampleInfo$sectors[[sector]]),
               list("decodedBy"=barcodesSample))
    }
    sampleInfo$callHistory <- append(sampleInfo$callHistory, match.call())
    decoded <- sampleInfo
    cleanit <- gc()
  } else {
    decoded <- splitByBarcode(sampleInfo, dnaSet, trimFrom=NULL, showStats=showStats, 
                              returnUnmatched=returnUnmatched, dereplicate=dereplicate)
    cleanit <- gc()
  }
  return(decoded)
}

#' @export
findBarcodes <- decodeByBarcode

#' Align a short pattern to variable length target sequences.
#'
#' Align a fixed length short pattern sequence (i.e. primers or adaptors) to subject sequences using \code{\link{pairwiseAlignment}}. This function uses default of type="overlap", gapOpening=-1, and gapExtension=-1 to align the patternSeq against subjectSeqs. One can adjust these parameters if prefered, but not recommended. This function is meant for aligning a short pattern onto large collection of subjects. If you are looking to align a vector sequence to subjects, then please use BLAT or see one of following \code{\link{blatSeqs}}, \code{\link{findAndRemoveVector}}
#'
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern. This is generally bigger than patternSeq, and cases where subjectSeqs is smaller than patternSeq will be ignored in the alignment.
#' @param patternSeq DNAString object or a sequence containing the query sequence to search. This is generally smaller than subjectSeqs. 
#' @param side which side of the sequence to perform the search: left, right or middle. Default is 'left'.
#' @param qualityThreshold percent of patternSeq to match. Default is 1, full match.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param bufferBases use x number of bases in addition to patternSeq length to perform the search. Beneficial in cases where the pattern has homopolymers or indels compared to the subject. Default is 5. Doesn't apply when side='middle'.
#' @param doRC perform reverse complement search of the defined pattern. Default is TRUE.
#' @param returnUnmatched return sequences which had no or less than 5\% match to the patternSeq. Default is FALSE.
#' @param returnLowScored return sequences which had quality score less than the defined qualityThreshold. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to FALSE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param ... extra parameters for \code{\link{pairwiseAlignment}}
#' @note 
#' \itemize{
#'  \item For qualityThreshold, the alignment score is calculated by (matches*2)-(mismatches+gaps) which programatically translates to round(nchar(patternSeq)*qualityThreshold)*2.
#'  \item Gaps and mismatches are weighed equally with value of -1 which can be overriden by defining extra parameters 'gapOpening' & 'gapExtension'.
#'  \item If qualityThreshold is 1, then it is a full match, if 0, then any match is accepted which is useful in searching linker sequences at 3' end. Beware, this function only searches for the pattern sequence in one orientation. If you are expecting to find the pattern in both orientation, you might be better off using BLAST/BLAT!
#'  \item If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#' }
#'
#' @return 
#' \itemize{
#'  \item IRanges object with starts, stops, and names of the aligned sequences. 
#'  \item If returnLowScored or returnUnmatched = T, then a CompressedIRangesList where x[["hits"]] has the good scoring hits, x[["Rejected"]] has the failed to match qualityThreshold hits, and x[["Absent"]] has the hits where the aligned bit is <=10\% match to the patternSeq.
#' }
#'
#' @seealso \code{\link{primerIDAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{doRCtest}}, \code{\link{findAndTrimSeq}}, \code{\link{blatSeqs}}, \code{\link{findAndRemoveVector}}
#'
#' @export
#'
#' @examples 
#' #pairwiseAlignSeqs(subjectSeqs, patternSeq, showStats=TRUE)
#' #pairwiseAlignSeqs(subjectSeqs, patternSeq, showStats=TRUE, qualityThreshold=0.5)
#'
pairwiseAlignSeqs <- function(subjectSeqs=NULL, patternSeq=NULL, side="left", 
                              qualityThreshold=1, showStats=FALSE, bufferBases=5, 
                              doRC=TRUE, returnUnmatched=FALSE, 
                              returnLowScored=FALSE, parallel=FALSE, ...) {
  
  .checkArgs_SEQed()
  
  ## only get the relevant side of subject sequence with extra bufferBases to 
  ## account for indels/mismatches & save memory while searching and avoid 
  ## searching elsewhere in the sequence
  if(tolower(side)=="left") {
    badSeqs <- DNAStringSet()
    culprits <- width(subjectSeqs) < (nchar(patternSeq)+bufferBases)
    if(any(culprits)) {
      badSeqs <- subjectSeqs[culprits]
      message(length(badSeqs),
              " sequences were removed from aligning since they were",
              " shorter than pattern getting aligned: ",
              (nchar(patternSeq)+bufferBases),"bp")            
      subjectSeqs <- subjectSeqs[!culprits]            
    }
    subjectSeqs2 <- subseq(subjectSeqs, start=1, end=(nchar(patternSeq)+bufferBases))
    overFromLeft <- rep(0,length(subjectSeqs))
  } else if (tolower(side)=="right") { 
    overFromLeft <- width(subjectSeqs)-(nchar(patternSeq)+bufferBases)
    overFromLeft[overFromLeft<1] <- 1
    subjectSeqs2 <- subseq(subjectSeqs, start=overFromLeft)
  } else {
    subjectSeqs2 <- subjectSeqs
    overFromLeft <- rep(0, length(subjectSeqs))
  }
  
  ## search both ways to test which side yields more hits!
  if(doRC) {
    patternSeq <- doRCtest(subjectSeqs2, patternSeq, qualityThreshold)            
  }
  
  if(parallel) {
    subjectSeqs2 <- chunkize(subjectSeqs2)
    hits <- foreach(x=iter(subjectSeqs2), .inorder=TRUE, 
                           .errorhandling="pass", 
                           .export=c("patternSeq",names(match.call())), 
                           .packages="Biostrings") %dopar% {
		if(any(names(match.call()) %in% c("type","gapOpening","gapExtension"))) {
		  pairwiseAlignment(x, patternSeq, ...)        
		} else {
		  pairwiseAlignment(x, patternSeq, type="overlap", 
							gapOpening=-1, gapExtension=-1, ...)
		}
    }
    hits <- do.call(c, hits)
  } else {
    ## type=overlap is best for primer trimming...see Biostrings Alignment vignette
    if(any(names(match.call()) %in% c("type","gapOpening","gapExtension"))) {
      hits <- pairwiseAlignment(subjectSeqs2, patternSeq, ...)        
    } else {
      hits <- pairwiseAlignment(subjectSeqs2, patternSeq, type="overlap", 
                                gapOpening=-1, gapExtension=-1, ...)
    }
  }
  
  stopifnot(length(hits)==length(subjectSeqs2))
  
  scores <- round(score(hits))
  highscored <- scores >= round(nchar(patternSeq)*qualityThreshold)*2
  unmatched <- nchar(hits) <= round(nchar(patternSeq)*.1) ## basically a small subset of highscored
  
  # no point in showing stats if all sequences are a potential match! #
  if(showStats & qualityThreshold!=0) {
    message("Total of ",as.numeric(table(highscored)['FALSE']),
            " did not have the defined pattern sequence (", patternSeq,
            ") that passed qualityThreshold on the ", side, " side")
  }
  
  ## extract starts-stops of the entire pattern hit ##
  starts <- start(pattern(hits))
  ends <- end(pattern(hits))
  namesq <- names(subjectSeqs)
  hits <- IRanges(start=starts+overFromLeft-ifelse(side=="right",2,0), 
                  end=ends+overFromLeft-ifelse(side=="right",2,0), 
                  names=namesq)
  rm("scores","subjectSeqs2","subjectSeqs","starts","ends","namesq")
  
  ## no need to test if there were any multiple hits since pairwiseAlignment will 
  ## only output one optimal alignment...see the man page.
  if(!returnLowScored & !returnUnmatched) {
    hits <- hits[highscored]
  } else {
    hitstoreturn <- IRangesList("hits"=hits[highscored])
    if(returnLowScored & length(hits[!highscored])>0) {
      hitstoreturn <- append(hitstoreturn, IRangesList("Rejected"=hits[!highscored]))
    }
    
    if(returnUnmatched & length(hits[unmatched])>0) {
      hitstoreturn <- append(hitstoreturn, IRangesList("Absent"=hits[unmatched]))
    }
    hits <- hitstoreturn
    rm(hitstoreturn)
  }             
  cleanit <- gc()
  return(hits)
}

#' Align a short pattern with PrimerID to variable length target sequences.
#'
#' Align a fixed length short pattern sequence containing primerID to variable length subject sequences using \code{\link{pairwiseAlignment}}. This function uses default of type="overlap", gapOpening=-1, and gapExtension=-1 to align the patterSeq against subjectSeqs. The search is broken up into as many pieces +1 as there are primerID and then compared against subjectSeqs. For example, patternSeq="AGCATCAGCANNNNNNNNNACGATCTACGCC" will launch two search jobs one per either side of Ns. For each search, qualityThreshold is used to filter out candidate alignments and the area in between is chosen to be the primerID. This strategy is benefical because of Indels introduced through homopolymer errors. Most likely the length of primerID(s) wont the same as you expected!
#'
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern. 
#' @param patternSeq DNAString object or a sequence containing the query sequence to search with the primerID.
#' @param qualityThreshold1 percent of first part of patternSeq to match. Default is 0.75.
#' @param qualityThreshold2 percent of second part of patternSeq to match. Default is 0.50.
#' @param doAnchored for primerID based patternSeq, use the base before and after primer ID in patternSeq as anchors?. Default is FALSE.
#' @param doRC perform reverse complement search of the defined pattern. Default is TRUE.
#' @param returnUnmatched return sequences if it had no or less than 5\% match to the first part of patternSeq before the primerID. Default is FALSE.
#' @param returnRejected return sequences if it only has a match to one side of patternSeq or primerID length does not match # of Ns +/-2 in the pattern. Default is FALSE.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param ... extra parameters for \code{\link{pairwiseAlignment}}
#'
#' @note 
#' \itemize{
#'  \item For qualityThreshold1 & qualityThreshold2, the alignment score is calculated by (matches*2)-(mismatches+gaps) which programatically translates to round(nchar(patternSeq)*qualityThreshold)*2
#'  \item Gaps and mismatches are weighed equally with value of -1 which can be overriden by defining extra parameters 'gapOpening' & 'gapExtension'.
#'  \item If qualityThreshold is 1, then it is a full match, if 0, then any match is accepted which is useful in searching linker sequences at 3' end. Beware, this function only searches for the pattern sequence in one orientation. If you are expecting to find the pattern in both orientation, you might be better off using BLAST/BLAT!
#' }
#'
#' @return 
#' \itemize{
#'  \item A CompressedIRangesList of length two, where x[["hits"]] is hits covering the entire patternSeq, and x[["primerIDs"]] is the potential primerID region. 
#'  \item If returnUnmatched = T, then x[["Absent"]] is appended which includes reads not matching the first part of patternSeq. 
#'  \item If returnRejected=TRUE, then x[["Rejected"]] includes reads that only matched one part of patternSeq or places where no primerID was found in between two part of patternSeq, and x[["RejectedprimerIDs"]] includes primerIDs that didn't match the correct length. 
#'  \item If doAnchored=TRUE, then x[["unAnchoredprimerIDs"]] includes reads that didn't match the base before and after primer ID on patternSeq.
#' }
#'
#' @seealso \code{\link{vpairwiseAlignSeqs}}, \code{\link{pairwiseAlignSeqs}}, \code{\link{doRCtest}}, \code{\link{blatSeqs}}, \code{\link{findAndRemoveVector}}
#'
#' @export
#'
#' @examples 
#' #primerIDAlignSeqs(subjectSeqs, patternSeq, showStats=TRUE)
#' #primerIDAlignSeqs(subjectSeqs, patternSeq, showStats=TRUE, qualityThreshold1=0.5)
#'
primerIDAlignSeqs <- function(subjectSeqs=NULL, patternSeq=NULL, 
                              qualityThreshold1=0.75, qualityThreshold2=0.50, 
                              doAnchored=FALSE, doRC=TRUE, returnUnmatched=FALSE, 
                              returnRejected=FALSE, showStats=FALSE, ...) {
  
  .checkArgs_SEQed()
  
  ## make sure there are Ns in the patternSeq for considering primerIDs
  if(length(unlist(gregexpr("N",patternSeq)))<4) {
    stop("There should be minimum of atleast 4 Ns in patternSeq to be ",
         "considered as a primerID sequence.")
  }
  
  ## Get the right orientation of the supplied patternSeq to peform proper search 
  ## at later two step search phase. 
  if(doRC) {
    patternSeq <- doRCtest(subjectSeqs, patternSeq)            
  }
  
  primerIDpos <- unlist(gregexpr("N", patternSeq))
  
  ## Perform primerID extraction by breaking the pattern into two parts...
  ## for sanity sakes due to homopolymers ##
  pattern1 <- as.character(subseq(DNAString(patternSeq), 1, primerIDpos[1]-1))
  pattern2 <- as.character(subseq(DNAString(patternSeq), primerIDpos[length(primerIDpos)]+1))
  
  pattern1.hits <- pairwiseAlignSeqs(subjectSeqs, pattern1, "middle", 
                                     qualityThreshold=qualityThreshold1, doRC=FALSE, 
                                     returnUnmatched=TRUE, ...)
  pattern2.hits <- pairwiseAlignSeqs(subjectSeqs, pattern2, "middle", 
                                     qualityThreshold=qualityThreshold2, doRC=FALSE, ...)
  
  ## Set aside reads which has no match to the pattern1...
  ## most likely mispriming if primerID is on 5' or read was too loong if on 3'
  ## Hits returned from pairwiseAlignSeqs will be filtered for low scored hits...
  ## so no need to check for those from subjectSeqs
  unmatched <- pattern1.hits[["Absent"]]
  pattern1.hits <- pattern1.hits[["hits"]]
  
  ## remove reads which only have a match to one of the patterns...crossover most likely!
  rejected1 <- setdiff(names(pattern1.hits), names(pattern2.hits))
  rejected2 <- setdiff(names(pattern2.hits), names(pattern1.hits))
  rejected <- c(pattern1.hits[names(pattern1.hits) %in% rejected1], 
                pattern2.hits[names(pattern2.hits) %in% rejected2])
  rm("rejected1","rejected2")
  if(showStats) { 
    message("Removed ", length(rejected),
            " read(s) for only matching one of pattern1 or pattern2") 
  }
  
  ## use only reads which match to both sides of the patterns.
  good.rows <- intersect(names(pattern1.hits), names(pattern2.hits))
  pattern1.hits <- pattern1.hits[names(pattern1.hits) %in% good.rows]
  pattern2.hits <- pattern2.hits[names(pattern2.hits) %in% good.rows]
  
  stopifnot(identical(names(pattern1.hits), names(pattern2.hits)))
  stopifnot(identical(names(pattern1.hits), good.rows))
  
  ## make sure there is no overlap of ranges between pattern1.hits & pattern2.hits
  ## if there is, then no primerID was found...remove it
  badAss <- end(pattern1.hits) >= start(pattern2.hits)
  if(any(badAss)) {
    rejected <- c(rejected, pattern1.hits[badAss], pattern2.hits[badAss])
    pattern1.hits <- pattern1.hits[!badAss]
    pattern2.hits <- pattern2.hits[!badAss]
    good.rows <- good.rows[!badAss]
    message("Removed ", table(badAss)["TRUE"],
            " read(s) for not having primerID present between pattern1 and pattern2")
  }
  
  hits <- IRanges(start=start(pattern1.hits), 
                  end=end(pattern2.hits), 
                  names=good.rows)        
  
  primerIDs <- IRanges(start=end(pattern1.hits)+1, 
                       end=start(pattern2.hits)-1, 
                       names=good.rows)        
  
  if(length(hits)==0) {
    stop("No hits found that matched both sides of patternSeq with ",
         "primerID in the middle.")
  }    
  
  ## do anchored search for only sequences that matched both sides of patternSeq
  if(doAnchored) {
    message("Found ", length(primerIDs), " total primerIDs before anchored filter.")
    
    ## get anchors of bases flanking Ns
    anchorBase.s <- substr(patternSeq, primerIDpos[1]-1, primerIDpos[1]-1)
    anchorBase.e <- substr(patternSeq, primerIDpos[length(primerIDpos)]+1, 
                           primerIDpos[length(primerIDpos)]+1)
    
    anchorBase.s.i <- trimSeqs(subjectSeqs, resize(pattern1.hits,width=1,fix="end"))
    anchorBase.e.i <- trimSeqs(subjectSeqs, resize(pattern2.hits,width=1,fix="start"))
    rows <- anchorBase.s==as.character(anchorBase.s.i) & anchorBase.e==as.character(anchorBase.e.i)
    
    unAnchored <- hits[!rows]
    primerIDs <- primerIDs[rows]
    hits <- hits[rows]
    message("Found ", length(primerIDs), " total primerIDs after anchored filter.")
    rm("rows","anchorBase.s.i","anchorBase.e.i")
  }
  rm("good.rows","pattern1.hits","pattern2.hits")
  cleanit <- gc()
  
  ## also remove any primerIDs that were too short or big than intended.
  nSize <- 2
  badAss <- !width(primerIDs) %in% (length(primerIDpos)-nSize):(length(primerIDpos)+nSize)
  if(any(badAss)) {
    rejectedprimerIDs <- hits[badAss]
    hits <- hits[!badAss]
    primerIDs <- primerIDs[!badAss]
    message("Removed ", table(badAss)["TRUE"],
            " read(s) for not having right primerID length")
  }
  
  hits <- IRangesList("hits"=hits, "primerIDs"=primerIDs)
  
  if(exists("unAnchored")) {
    if(length(unAnchored)>0) { 
      hits <- append(hits, IRangesList("unAnchoredprimerIDs"=unAnchored)) 
    }
  }
  
  if(returnUnmatched & length(unmatched)>0) {
    hits <- append(hits, IRangesList("Absent"=unmatched))
  }
  
  if(returnRejected) {
    if(length(rejected)>0) { 
      hits <- append(hits, IRangesList("Rejected"=rejected)) 
    }
    if(exists("rejectedprimerIDs")) { 
      if(length(rejectedprimerIDs)>0) { 
        hits <- append(hits, IRangesList("RejectedprimerIDs"=rejectedprimerIDs)) 
      } 
    }
  }    
  return(hits)
}

#' Align a short pattern to variable length target sequences.
#'
#' Align a fixed length short pattern sequence to subject sequences using \code{\link{vmatchPattern}}. This function is meant for aligning a short pattern onto large collection of subjects. If you are looking to align a vector sequence to subjects, then please use BLAT.
#'
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern. This is generally bigger than patternSeq, and cases where subjectSeqs is smaller than patternSeq will be ignored in the alignment.
#' @param patternSeq DNAString object or a sequence containing the query sequence to search. This is generally smaller than subjectSeqs. 
#' @param side which side of the sequence to perform the search: left, right, or middle. Default is 'left'.
#' @param qualityThreshold percent of patternSeq to match. Default is 1, full match. This is supplied to max.mismatch parameter of \code{\link{vmatchPattern}} as round(nchar(patternSeq)*(1-qualityThreshold)).
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param bufferBases use x number of bases in addition to patternSeq length to perform the search. Beneficial in cases where the pattern has homopolymers or indels compared to the subject. Default is 5. Doesn't apply when side='middle'.
#' @param doRC perform reverse complement search of the defined pattern. Default is TRUE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to FALSE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param ... extra parameters for \code{\link{vmatchPattern}} except for 'max.mismatch' since it's calculated internally using the 'qualityThreshold' parameter.
#'
#'
#' @note 
#' \itemize{
#'  \item For qualityThreshold, the alignment score is calculated by (matches*2)-(mismatches+gaps) which programatically translates to round(nchar(patternSeq)*qualityThreshold)*2.
#'  \item No indels are allowed in the function, if expecting indels then use \code{\link{pairwiseAlignSeqs}}.
#'  \item If qualityThreshold is 1, then it is a full match, if 0, then any match is accepted which is useful in searching linker sequences at 3' end. Beware, this function only searches for the pattern sequence in one orientation. If you are expecting to find the pattern in both orientation, you might be better off using BLAST/BLAT!
#'  \item If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#' }
#' 
#' @return IRanges object with starts, stops, and names of the aligned sequences.
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{doRCtest}}, \code{\link{findAndTrimSeq}}, \code{\link{blatSeqs}}, \code{\link{findAndRemoveVector}}
#'
#' @export
#'
#' @examples 
#' #vpairwiseAlignSeqs(subjectSeqs, patternSeq, showStats=TRUE)
#' #vpairwiseAlignSeqs(subjectSeqs, patternSeq, showStats=TRUE, qualityThreshold=0.5)
#'
vpairwiseAlignSeqs <- function(subjectSeqs=NULL, patternSeq=NULL, side="left", 
                               qualityThreshold=1, showStats=FALSE, 
                               bufferBases=5, doRC=TRUE, parallel=FALSE, ...) {
  
  .checkArgs_SEQed()
  
  ## Only get the relevant side of subject sequence with extra bufferBases to...
  ## account for indels/mismatches & save memory while searching & avoid searching...
  ## elsewhere in the sequence
  if(tolower(side)=="left") {
    badSeqs <- DNAStringSet()
    culprits <- width(subjectSeqs) < (nchar(patternSeq)+bufferBases)
    if(any(culprits)) {
      badSeqs <- subjectSeqs[culprits]
      message(length(badSeqs),
              " sequences were removed from aligning since they were",
              " shorter than pattern getting aligned: ",
              (nchar(patternSeq)+bufferBases),"bp")
      subjectSeqs <- subjectSeqs[!culprits]            
    }
    subjectSeqs2 <- subseq(subjectSeqs, start=1, end=(nchar(patternSeq)+bufferBases))
    overFromLeft <- rep(0, length(subjectSeqs))
  } else if (tolower(side)=="right") { 
    overFromLeft <- width(subjectSeqs)-(nchar(patternSeq)+bufferBases)
    overFromLeft[overFromLeft<1] <- 1
    subjectSeqs2 <- subseq(subjectSeqs, start=overFromLeft)
  } else {
    subjectSeqs2 <- subjectSeqs
    overFromLeft <- rep(0, length(subjectSeqs))
  }
  
  ## search both ways to test which side yields more hits!        
  if(doRC) {        
    patternSeq <- doRCtest(subjectSeqs2, patternSeq, qualityThreshold)            
  }
  
  if(parallel) {
    subjectSeqs2 <- chunkize(subjectSeqs2)
    hits <- foreach(x=iter(subjectSeqs2), .inorder=TRUE, 
                    .errorhandling="pass", 
                    .export=c("patternSeq", names(match.call())), 
                    .packages="Biostrings") %dopar% {
                      maxmis <- round(nchar(patternSeq)*(1-qualityThreshold))
                      bore <- vmatchPattern(patternSeq, x, 
                                            max.mismatch=maxmis, ...)
                      unlist(bore, recursive=TRUE, use.names=TRUE)              
                    }
    hits <- do.call(c, hits)
  } else {
    hits <- vmatchPattern(patternSeq, subjectSeqs2, 
                        max.mismatch=round(nchar(patternSeq)*(1-qualityThreshold)), ...)
  	hits <- unlist(hits, recursive=TRUE, use.names=TRUE)
  }
    
  ## test if there were any multiple hits which are overlapping and if they are reduce them
  counts <- Rle(names(hits))
  if(any(runLength(counts)>1)) {
    reduced <- reduce(GRanges(seqnames=names(hits),
                              IRanges(start=start(hits), end=end(hits))))
    counts <- seqnames(reduced)
    hits <- ranges(reduced)
    names(hits) <- as.character(seqnames(reduced))
    rm(reduced)
    if(any(runLength(counts)>1)) {
      message("More than 1 pattern (", patternSeq,") match found for:",
              paste(runValue(counts)[runLength(counts)>1],collapse=","),
              "\nUsing the latter occuring hit as the dominant for each read.")
      toremove <- c()
      for(culprits in as.character(runValue(counts)[runLength(counts)>1])) {
        rows <- which(names(hits) %in% culprits)
        toremove <- c(toremove, rows[1:length(rows)-1])
      }
      hits <- hits[-toremove]
      counts <- Rle(names(hits))
      if(any(runLength(counts)>1)) {
        stop("More than 1 pattern unresolved (",patternSeq,") match found for:",
             paste(runValue(counts)[runLength(counts)>1],collapse=","))
      }
    }
  }
  
  good.row <- names(subjectSeqs2) %in% names(hits)
  
  if(showStats) {
    message("Total of ", as.numeric(table(good.row)['FALSE']),
            " did not have the defined pattern sequence (", patternSeq,
            ") that passed qualityThreshold on the ", side, " side")
  }  
  
  starts <- start(hits)
  ends <- end(hits)
  namesq <- names(hits)
  rm("hits","subjectSeqs2")     
  
  hits <- IRanges(start=starts+overFromLeft[good.row]-ifelse(side=="right",2,0), 
                  end=ends+overFromLeft[good.row]-ifelse(side=="right",2,0), 
                  names=namesq)
  
  cleanit <- gc()
  return(hits)
}

#' Test if pattern aligns better in +/- orientation.
#'
#' Given a fixed length pattern sequence and variable length subject sequences, the function roughly finds which orientation of pattern yields the most hits. The function doing the heavylifting is \code{\link{vcountPattern}}. This is an accessory function used in function listed under See Also section below. 
#'
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern. 
#' @param patternSeq DNAString object or a sequence containing the query sequence to search.
#' @param qualityThreshold percent of patternSeq to match. Default is 0.50, half match. This is supplied to max.mismatch parameter of \code{\link{vcountPattern}} as round(nchar(patternSeq)*(1-qualityThreshold)).
#' @param cores.use Number of cores to use for aligning. Only need two max! Default is 2. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @return patternSeq that aligned the best. 
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{primerIDAlignSeqs}}
#'
#' @export
#'
#' @examples 
#' #doRCtest(subjectSeqs,patternSeq)
#'
doRCtest <- function(subjectSeqs=NULL, patternSeq=NULL, 
                     qualityThreshold=0.5, core.use=2) {
  
  .checkArgs_SEQed()
  
  if(core.use==1) { registerDoSEQ() }
  
  patternSeq.rc <- as.character(reverseComplement(DNAString(patternSeq)))
  hits <- foreach(x=iter(c(patternSeq,patternSeq.rc)), .inorder=TRUE, 
                  .errorhandling="pass", 
                  .export=c("subjectSeqs","qualityThreshold"), 
                  .packages="Biostrings") %dopar% {
    counts <- pmin(vcountPattern(x, subjectSeqs, 
                                 max.mismatch=round(nchar(x)*(1-qualityThreshold))),1)
    sum(counts)
  }
  
  if(all(hits==0)) {
    stop("No hits found")
  }
  
  if(hits[[1]] < hits[[2]]) {
    message("There were less/no good hits found for pattern ", patternSeq,
            "than its reverse complement...", patternSeq.rc,
            "\nUsing the latter to perform the searching.")
    patternSeq <- patternSeq.rc
  }
  
  return(patternSeq)
}

#' This function generates alignment statistics for findXXXXXX functions. Evaluation of this function happens in the parent function.
#'
.showFindStats <- function() {
  checks <- expression(
    isFindLinker <- grepl("linker",trimmedObj,ignore.case=TRUE) |
                    grepl("linker",featureTrim,ignore.case=TRUE),
    listed <- sapply(get(trimmedObj), is.list),
    if(any(listed)) {
      totals <- if(isFindLinker) {        
        sapply(get(trimmedObj)[listed], sapply, function(x) length(x$hits))
      } else {
        sapply(get(trimmedObj)[listed], sapply, length)
      }
      counts <- sapply(get(rawObj)[names(get(trimmedObj)[listed])], sapply, length)
      stopifnot(identical(colnames(totals), colnames(counts)))
      toprint <- cbind(data.frame("Total"=t(totals)),
                       data.frame("valueTempCol"=t(100*(totals/counts))))
      names(toprint)[3:4] <- gsub("valueTempCol", valueColname, names(toprint)[3:4])
      mean.test <- mean(toprint[,paste0(valueColname,
                                        ifelse(isFindLinker,".pair2",".pair1"))])<=5      
    } else {
      if(isFindLinker) {
        toprint <- data.frame("Total"=sapply(sapply(linkerTrimmed[!listed], 
                                                    "[[", "hits"), length))
      } else {
        toprint <- as.data.frame(sapply(get(trimmedObj)[!listed], length))
      }
      names(toprint) <- "Total"            
      counts <- sapply(get(rawObj)[names(get(trimmedObj)[!listed])], length)
      toprint[,valueColname] <- 100*(toprint$Total/counts[rownames(toprint)])
      mean.test <- mean(toprint[,valueColname])<=5
    },
    
    toprint$SampleName <- rownames(toprint),
    rownames(toprint) <- NULL,
    
    if(showStats) {            
      message("Sequence Stats after ",featureTrim," alignment:")
      print(toprint)        
    },
    
    ## if <= 5% of sequences found primers...then something is wrong with the primer sequences provided???
    if(mean.test) {
      stop("Something seems to be wrong with the ",featureTrim,"s provided for ",
           "each sample. On average <= 5% of sequences found ",featureTrim," match ",
           "for the entire run!!!")
    }
  )
  
  eval.parent(checks)
}

#' Find the 5' primers and add results to SampleInfo object. 
#'
#' Given a sampleInfo object, the function finds 5' primers for each sample per sector and adds the results back to the object. This is a specialized function which depends on many other functions shown in 'see also section' to perform specialized trimming of 5' primer/adaptor found in the sampleInfo object. The sequence itself is never trimmed but rather coordinates of primer portion is recorded back to the object and used subsequently by \code{\link{extractSeqs}} function to perform the trimming.
#'
#' @param sampleInfo sample information SimpleList object outputted from \code{\link{decodeByBarcode}}, which holds decoded sequences for samples per sector/quadrant along with information of sample to primer associations.
#' @param alignWay method to utilize for detecting the primers. One of following: "slow" (Default), or "fast". Fast, calls \code{\link{vpairwiseAlignSeqs}} and uses \code{\link{vpatternMatch}} at its core, which is less accurate with indels and mismatches but much faster. Slow, calls \code{\link{pairwiseAlignSeqs}} and uses \code{\link{pairwiseAlignment}} at its core, which is accurate with indels and mismatches but slower.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param doRC perform reverse complement search of the defined pattern/primer. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param samplenames a vector of samplenames to process. Default is NULL, which processes all samples from sampleInfo object.
#' @param bypassChecks skip checkpoints which detect if something was odd with the data? Default is FALSE.
#' @param ... extra parameters to be passed to either \code{\link{vpatternMatch}} or \code{\link{pairwiseAlignment}} depending on 'alignWay' parameter.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new data added under each sector and sample. New data attributes include: primed
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{extractFeature}}, \code{\link{extractSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}}, \code{\link{findLinkers}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #findPrimers(sampleInfo, showStats=TRUE)
#' #findPrimers(sampleInfo, alignWay="slow", showStats=TRUE)
#'
findPrimers <- function(sampleInfo, alignWay="slow", showStats=FALSE, 
                        doRC=FALSE, parallel=TRUE, samplenames=NULL, 
                        bypassChecks=FALSE, ...) {    
  
  .checkArgs_SEQed()
  
  ## test if there are decoded sequences in the sampleinfo object ##
  decoded <- extractFeature(sampleInfo, feature="decoded")
  samplesDecoded <- sapply(decoded, names, simplify=FALSE)
  sectorsDecoded <- names(which(sapply(sapply(decoded,length), ">", 0)))
  rm(decoded)
  cleanit <- gc()
  
  if(length(sectorsDecoded)==0) {
    stop("No decoded information found in sampleInfo...",
         "did you run decodeByBarcode()?")
  }
  
  for(sector in sectorsDecoded) {
    message("Processing sector ", sector)
    
    ## prepare sample to primer association ##
    ltrPrimers <- toupper(extractFeature(sampleInfo, sector=sector,
                                         feature="primerltrsequence")[[sector]])
    skippers <- ltrPrimers=="SKIP"
    if(!all(skippers) & (length(ltrPrimers[!skippers])==0 | 
                           mean(nchar(ltrPrimers[!skippers]))<=5 | 
                           all(is.na(ltrPrimers[!skippers])))) {
      stop("Either the primer size is too short (<=5) or ",
           "no primers are found in sample information object.")
    }
    
    ## refine sample list if specific samples are supplied ##
    samplesToProcess <- samplesDecoded[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }
    
    ## find any samples which need to be skipped
    skip.samples <- names(which(skippers))
    if(length(skip.samples)>0) {
      samplesToProcess <- samplesToProcess[!samplesToProcess %in% skip.samples]
      message("Skipping samples ", paste(skip.samples,collapse=","))
      sampleInfo <- addFeature(sampleInfo, sector, skip.samples, feature="primed", 
                               value=structure(rep("SKIPPED",length(skip.samples)), 
                                               names=skip.samples))
    }

    ## dont bother searching if no samples are to be processed! ##
    if(length(samplesToProcess)>0) {			
      ## get the decoded reads ##
      decoded <- extractSeqs(sampleInfo, sector, samplesToProcess,
                             feature="decoded")[[sector]]
      
      ## find paired end samples
      isPaired <- extractFeature(sampleInfo, sector, feature='pairedend')[[sector]]
      
      ## trim the primers ##
      message("\tFinding Primers.")
      primerIdentity <- extractFeature(sampleInfo, sector,
                                       feature="primerltridentity")[[sector]]
      stopifnot(length(primerIdentity)>0)
      
      primerTrimmed <- foreach(x=iter(samplesToProcess), .inorder=TRUE, 
                               .errorhandling="pass", 
                               .export=c("ltrPrimers", "primerIdentity", "doRC",
                                         "alignWay", "vpairwiseAlignSeqs", 
                                         "pairwiseAlignSeqs", "decoded", "isPaired"), 
                               .packages="Biostrings") %dopar% {
        if(isPaired[[x]]) {
          p1 <- switch(alignWay,
                       fast = vpairwiseAlignSeqs(decoded[[x]]$pair1, 
                                                 ltrPrimers[[x]], "left", 
                                                 (primerIdentity[[x]]-.05), 
                                                 doRC=doRC, ...),
                       slow = pairwiseAlignSeqs(decoded[[x]]$pair1, 
                                                ltrPrimers[[x]], "left", 
                                                (primerIdentity[[x]]), 
                                                doRC=doRC, ...)                
          )
          
          ## side needs to be middle since we dont know where in pair2 the primer can be
          p2 <- switch(alignWay,
                       fast = vpairwiseAlignSeqs(decoded[[x]]$pair2, 
                                                 ltrPrimers[[x]], "middle", 
                                                 (primerIdentity[[x]]-.05), 
                                                 doRC=doRC, ...),
                       slow = pairwiseAlignSeqs(decoded[[x]]$pair2, 
                                                ltrPrimers[[x]], "middle", 
                                                (primerIdentity[[x]]), 
                                                doRC=doRC, ...)                
          )
          list("pair1"=p1, "pair2"=p2)
        } else {
          switch(alignWay,
                 fast = vpairwiseAlignSeqs(decoded[[x]], ltrPrimers[[x]], "left", 
                                           qualityThreshold=(primerIdentity[[x]]-.05), 
                                           doRC=doRC, ...),
                 slow = pairwiseAlignSeqs(decoded[[x]], ltrPrimers[[x]], "left", 
                                          qualityThreshold=(primerIdentity[[x]]), 
                                          doRC=doRC, ...)                
          )
        }
      }
      names(primerTrimmed) <- samplesToProcess
      
      ## check if any error occured during alignments ##
      if(any(grepl("simpleError",primerTrimmed))) {
        stop("Error encountered in LTR Trimming function",
             paste(names(primerTrimmed[grepl("simpleError",primerTrimmed)]),
                   collapse=", "))
      }        
      
      ## remove samples with no primer hits from further processing ##
      culprits <- grep("No hits found",primerTrimmed)
      if(length(culprits)>0) {
        message("Following sample(s) had no hits for primer alignment: ",
                paste(samplesToProcess[culprits],collapse=", "))
        samplesToProcess <- samplesToProcess[-c(culprits)]
        primerTrimmed <- primerTrimmed[-c(culprits)]
      }
      
      cleanit <- gc()
      
      if(!bypassChecks | showStats) {
        eval(expression(trimmedObj <- "primerTrimmed", rawObj <- "decoded", 
                        featureTrim <- "Primer", valueColname <- "PercOfDecoded"))
        .showFindStats()
      }
      
      ## modify metadata attribute, write primer coordinates back to sampleInfo object & trim
      message("Adding primer info back to the object")
      sampleInfo <- addFeature(sampleInfo, sector, names(primerTrimmed),
                               feature="primed", value=primerTrimmed)
      rm(primerTrimmed)
      cleanit <- gc()
    }
  }
  sampleInfo$callHistory <- append(sampleInfo$callHistory, match.call())
  return(sampleInfo)
}

#' Find the 5' LTRs and add results to SampleInfo object. 
#'
#' Given a sampleInfo object, the function finds 5' LTR following the primer for each sample per sector and adds the results back to the object. This is a specialized function which depends on many other functions shown in 'see also section' to perform specialized trimming of 5' viral LTRs found in the sampleInfo object. The sequence itself is never trimmed but rather coordinates of LTR portion is added to primer coordinates and recorded back to the object and used subsequently by \code{\link{extractSeqs}} function to perform the trimming. This function heavily relies on \code{\link{pairwiseAlignSeqs}}.
#'
#' @param sampleInfo sample information SimpleList object outputted from \code{\link{findPrimers}}, which holds decoded and primed sequences for samples per sector/quadrant along with information of sample to LTR associations.
#' @param showStats toggle output of search statistics. Default is FALSE. For paired end data, stats for "pair2" is relative to decoded and/or primed reads.
#' @param doRC perform reverse complement search of the defined pattern/LTR sequence. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param samplenames a vector of samplenames to process. Default is NULL, which processes all samples from sampleInfo object.
#' @param bypassChecks skip checkpoints which detect if something was odd with the data? Default is FALSE.
#' @param ... extra parameters to be passed to \code{\link{pairwiseAlignment}}.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new data added under each sector and sample. New data attributes include: LTRed
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{extractFeature}}, \code{\link{extractSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findPrimers}}, \code{\link{findLinkers}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #findLTRs(sampleInfo,showStats=TRUE)
#'
findLTRs <- function(sampleInfo, showStats=FALSE, doRC=FALSE, 
                     parallel=TRUE, samplenames=NULL, bypassChecks=FALSE, ...) {    
  
  .checkArgs_SEQed()
  
  ## test if there are primed sequences in the sampleinfo object ##   
  primed <- extractFeature(sampleInfo, feature="primed")
  samplesprimed <- sapply(primed, names, simplify=FALSE)
  sectorsPrimed <- names(which(sapply(sapply(primed, length), ">", 0)))
  rm(primed)
  cleanit <- gc()
  
  if(length(sectorsPrimed)==0) {
    stop("No primed information found in sampleInfo...did you run findPrimers()?")
  }
  
  for(sector in sectorsPrimed) {
    message("Processing sector ",sector)
    
    ## prepare sample to LTR bit associations ##
    sampleLTRbits <- toupper(extractFeature(sampleInfo,sector=sector,
                                            feature="ltrbitsequence")[[sector]])
    skippers <- sampleLTRbits=="SKIP"
    if(!all(skippers) & (length(sampleLTRbits[!skippers])==0 | 
                           mean(nchar(sampleLTRbits[!skippers]))<=1 | 
                           all(is.na(sampleLTRbits[!skippers])))) {
      stop("Either LTR bit sequence is too short (<=1) or no LTR bits ",
           "found in sample information object.")
    }
    
    ## refine sample list if specific samples are supplied ##
    samplesToProcess <- samplesprimed[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }
    
    # find any samples which need to be skipped
    skip.samples <- names(which(skippers))
    if(length(skip.samples)>0) {
      samplesToProcess <- samplesToProcess[!samplesToProcess %in% skip.samples]
      message("Skipping samples ", paste(skip.samples,collapse=","))
      sampleInfo <- addFeature(sampleInfo, sector, skip.samples, 
                               feature="LTRed", 
                               value=structure(rep("SKIPPED",length(skip.samples)), 
                                               names=skip.samples))
    }
    
    ## dont bother searching if no samples are to be processed! ##
    if(length(samplesToProcess)>0) {
      ## get the primer trimmed reads ##
      primerTrimmed <- extractSeqs(sampleInfo, sector, samplesToProcess,
                                   feature="primed")[[sector]]
      
      ## find paired end samples...
      ## add reads which aren't primed in pair2 for cases where reads were long, etc.
      isPaired <- extractFeature(sampleInfo, sector, samplesToProcess,
                                 feature='pairedend')[[sector]]
      if(any(isPaired)) {
        message("Getting reads from pair2 which weren't primed...")
        decoded <- extractSeqs(sampleInfo, sector, names(which(isPaired)),
                               feature="decoded", pairReturn='pair2')[[sector]]
        rows <- intersect(names(decoded), names(primerTrimmed))
        
        bore <- mapply(function(x,y) {
          x$pair2 <- c(x$pair2, y[!names(y) %in% names(x$pair2)])
          x
        }, primerTrimmed[rows], decoded[rows])
        primerTrimmed[rows] <- as(bore, "DataFrame")        
        rm(bore)
      }
      
      ## trim LTRbits using slow method since its the best! ##
      message("\tFinding LTR bits.")
      ltrBitIdentity <- extractFeature(sampleInfo,sector=sector,
                                       feature="ltrbitidentity")[[sector]]
      
      ltrTrimmed <- foreach(x=iter(samplesToProcess), .inorder=TRUE, 
                            .errorhandling="pass", 
                            .export=c("primerTrimmed", "sampleLTRbits", "isPaired",
                                      "ltrBitIdentity", "doRC", "pairwiseAlignSeqs"), 
                            .packages="Biostrings") %dopar% {
        if(isPaired[[x]]) {
          p1 <- pairwiseAlignSeqs(primerTrimmed[[x]]$pair1, sampleLTRbits[[x]], "left", 
                                  qualityThreshold=ltrBitIdentity[[x]], doRC=doRC, ...)
          
          p2 <- pairwiseAlignSeqs(primerTrimmed[[x]]$pair2, sampleLTRbits[[x]], 
                                  "middle", 
                                  qualityThreshold=ltrBitIdentity[[x]], doRC=doRC, ...)

          list("pair1"=p1, "pair2"=p2)          
        } else {
          pairwiseAlignSeqs(primerTrimmed[[x]], sampleLTRbits[[x]], "left", 
                            qualityThreshold=ltrBitIdentity[[x]], doRC=doRC, ...)
        }        
      }
      names(ltrTrimmed) <- samplesToProcess
      
      ## check if any error occured during alignments ##
      if(any(grepl("simpleError",ltrTrimmed))) {
        stop("Error encountered in LTR Trimming function",
             paste(names(ltrTrimmed[grepl("simpleError",ltrTrimmed)]),
                   collapse=", "))
      }
      
      ## remove samples with no LTR hits from further processing ##
      culprits <- grep("No hits found",ltrTrimmed)
      if(length(culprits)>0) {
        message("Following sample(s) had no hits for LTR bit alignment: ",
                paste(samplesToProcess[culprits],collapse=", "))
        samplesToProcess <- samplesToProcess[-c(culprits)]
        ltrTrimmed <- ltrTrimmed[-c(culprits)]
      }    

      cleanit <- gc()
      
      if(!bypassChecks | showStats) {
        eval(expression(trimmedObj <- "ltrTrimmed", rawObj <- "primerTrimmed", 
                        featureTrim <- "LTR", valueColname <- "PercOfPrimed"))
        .showFindStats()
      }
      
      ## modify metadata attribute, add LTR bit coordinates to primer and write back to sampleInfo object & trim...remember everything is relative to the entire read length!
      message("Adding LTR info back to the object")
      for(x in names(ltrTrimmed)) {
        cat(".")
        if(isPaired[[x]]) {
          worked <- sapply(ltrTrimmed[[x]],length)>0
          if(any(worked)) {
            for(pair in names(which(worked))) {
              primed <- sampleInfo$sectors[[sector]]$samples[[x]]$primed[[pair]]
              rows <- names(ltrTrimmed[[x]][[pair]]) %in% names(primed)
              primed.end <- end(primed[names(ltrTrimmed[[x]][[pair]][rows])])
              rm(primed)
              
              end(ltrTrimmed[[x]][[pair]][rows]) <- 
                end(ltrTrimmed[[x]][[pair]][rows]) + primed.end              
              start(ltrTrimmed[[x]][[pair]][rows]) <- 
                start(ltrTrimmed[[x]][[pair]][rows]) + primed.end              
              rm(primed.end)
            }
          }
          sampleInfo$sectors[[sector]]$samples[[x]]$LTRed <- ltrTrimmed[[x]]
        } else {
          if(length(ltrTrimmed[[x]])>0) {
            primed <- sampleInfo$sectors[[sector]]$samples[[x]]$primed
            primed.end <- end(primed[names(ltrTrimmed[[x]])])
            rm(primed)
            
            end(ltrTrimmed[[x]]) <- end(ltrTrimmed[[x]]) + primed.end
            start(ltrTrimmed[[x]]) <- start(ltrTrimmed[[x]]) + primed.end
            rm(primed.end)
            
            sampleInfo$sectors[[sector]]$samples[[x]]$LTRed <- ltrTrimmed[[x]]
          }
        }        
      }
      rm("ltrTrimmed","primerTrimmed")
      cleanit <- gc()
    }
  }
  sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
  return(sampleInfo)
}

#' Find the 3' linkers and add results to SampleInfo object. 
#'
#' Given a sampleInfo object, the function finds 3' linkers for each sample per sector and adds the results back to the object. This is a specialized function which depends on many other functions shown in 'see also section' to perform specialized trimming of 3' primer/linker adaptor sequence found in the sampleInfo object. The sequence itself is never trimmed but rather coordinates of linker portion is recorded back to the object and used subsequently by \code{\link{extractSeqs}} function to perform the trimming. This function heavily relies on either \code{\link{pairwiseAlignSeqs}} or \code{\link{primerIDAlignSeqs}} depending upon whether linkers getting aligned have primerID in it or not.
#'
#' @param sampleInfo sample information SimpleList object outputted from \code{\link{findPrimers}} or \code{\link{findLTRs}}, which holds decoded sequences for samples per sector/quadrant along with information of sample to primer associations.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param doRC perform reverse complement search of the defined pattern/linker sequence. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param samplenames a vector of samplenames to process. Default is NULL, which processes all samples from sampleInfo object.
#' @param bypassChecks skip checkpoints which detect if something was odd with the data? Default is FALSE.
#' @param ... extra parameters to be passed to \code{\link{pairwiseAlignment}}.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new data added under each sector and sample. New data attributes include: linkered. If linkers have primerID then, primerIDs attribute is appended as well. 
#'
#' @note If no linker matches are found with default options, then try doRC=TRUE. If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}}, \code{\link{findPrimers}}, \code{\link{extractFeature}}, \code{\link{extractSeqs}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #findLinkers(sampleInfo,showStats=TRUE)
#'
findLinkers <- function(sampleInfo, showStats=FALSE, doRC=FALSE, parallel=TRUE, 
                        samplenames=NULL, bypassChecks=FALSE, ...) {    
  
  .checkArgs_SEQed()
  
  ## test if there are decoded sequences in the sampleinfo object ##
  decoded <- extractFeature(sampleInfo, feature="decoded")
  toProcessSamples <- sapply(decoded, names, simplify=FALSE)
  sectorsDecoded <- names(which(sapply(sapply(decoded,length), ">", 0)))
  rm(decoded)
  cleanit <- gc()
  
  if(length(sectorsDecoded)==0) {
    stop("No decoded information found in sampleInfo...",
         "did you run findBarcodes()/decodeByBarcode()?")
  }
  
  for(sector in sectorsDecoded) {
    message("Processing sector ",sector)
    
    ## prepare sample to linker association ##
    sampleLinkers <- toupper(extractFeature(sampleInfo, sector=sector,
                                            feature="linkersequence")[[sector]])
    skippers <- sampleLinkers=="SKIP"
    if(!all(skippers) & (length(sampleLinkers[!skippers])==0 | 
                           mean(nchar(sampleLinkers[!skippers]))<=10 | 
                           all(is.na(sampleLinkers[!skippers])))) {
      stop("Either Linker sequence is too short (<=10) or ",
           "no Linkers found in sample information object.")
    }
    
    ## get the linker quality thresholds for non primerID based samples ##
    linkerIdentity <- extractFeature(sampleInfo, sector=sector,
                                     feature="linkeridentity")[[sector]]
    stopifnot(length(linkerIdentity)>0)
    
    ## check if any are primerIDed and get their identity thresholds ##
    primerIded <- extractFeature(sampleInfo, sector=sector,
                                 feature="primeridinlinker")[[sector]]        
    primerIded.threshold1 <- extractFeature(sampleInfo, sector=sector,
                                            feature="primeridinlinkeridentity1")[[sector]]        
    primerIded.threshold2 <- extractFeature(sampleInfo, sector=sector, 
                                            feature="primeridinlinkeridentity2")[[sector]]        
    
    ## subset samplenames from all samples if defined ##
    samplesToProcess <- toProcessSamples[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }
    
    # find any samples which need to be skipped
    skip.samples <- names(which(skippers))
    if(length(skip.samples)>0) {
      samplesToProcess <- samplesToProcess[!samplesToProcess %in% skip.samples]
      message("Skipping samples ", paste(skip.samples,collapse=","))
      sampleInfo <- addFeature(sampleInfo, sector, skip.samples, feature="linkered", 
                               value=structure(rep("SKIPPED",length(skip.samples)), 
                                               names=skip.samples))
    }
    
    ## dont bother searching if no samples are to be processed! ##
    if(length(samplesToProcess)>0) {
      
      toProcess <- extractSeqs(sampleInfo, sector, samplesToProcess,
                               feature="decoded")[[sector]]
      
      ## find paired end samples
      isPaired <- extractFeature(sampleInfo, sector, feature='pairedend')[[sector]]
      
      ## trim the Linkers ##
      message("\tFinding Linkers.")
      linkerTrimmed <- foreach(x=iter(samplesToProcess), .inorder=TRUE, 
                               .errorhandling="pass", 
                               .export=c("primerIded", "toProcess", "sampleLinkers",
                                         "doRC", "linkerIdentity", "isPaired",
                                         "primerIded.threshold1",
                                         "primerIded.threshold2",
                                         "pairwiseAlignSeqs",
                                         "primerIDAlignSeqs"), 
                               .packages="Biostrings") %dopar% {        
        if(primerIded[[x]]) {
          if(isPaired[[x]]) {
            p1 <- primerIDAlignSeqs(toProcess[[x]]$pair1, sampleLinkers[[x]], 
                                    doAnchored=TRUE, returnUnmatched=TRUE, 
                                    returnRejected=TRUE, doRC=doRC, 
                                    qualityThreshold1=primerIded.threshold1[[x]], 
                                    qualityThreshold2=primerIded.threshold2[[x]], ...)
            
            p2 <- primerIDAlignSeqs(toProcess[[x]]$pair2, sampleLinkers[[x]], 
                                    doAnchored=TRUE, returnUnmatched=TRUE, 
                                    returnRejected=TRUE, doRC=doRC, 
                                    qualityThreshold1=primerIded.threshold1[[x]], 
                                    qualityThreshold2=primerIded.threshold2[[x]], ...)
            
            list("pair1"=p1, "pair2"=p2)
          } else {
            primerIDAlignSeqs(toProcess[[x]], sampleLinkers[[x]], doAnchored=TRUE, 
                              returnUnmatched=TRUE, returnRejected=TRUE, doRC=doRC, 
                              qualityThreshold1=primerIded.threshold1[[x]], 
                              qualityThreshold2=primerIded.threshold2[[x]], ...)
          }          
        } else {
          ## use side="middle" since more junk sequence can be present after linker which would fail pairwiseAlignSeqs if side='right' for single end reads or "pair1"
          if(isPaired[[x]]) {
            p1 <- pairwiseAlignSeqs(toProcess[[x]]$pair1, sampleLinkers[[x]], "middle", 
                                    qualityThreshold=linkerIdentity[[x]], 
                                    returnUnmatched=TRUE, returnLowScored=TRUE, 
                                    doRC=doRC, ...)
            
            # pair2 should end with linker, hence side='right'            
            p2 <- pairwiseAlignSeqs(toProcess[[x]]$pair2, sampleLinkers[[x]], "right", 
                                    qualityThreshold=linkerIdentity[[x]], 
                                    returnUnmatched=TRUE, returnLowScored=TRUE, 
                                    doRC=doRC, ...)
            
            list("pair1"=p1, "pair2"=p2)
          } else {
            pairwiseAlignSeqs(toProcess[[x]], sampleLinkers[[x]], "middle", 
                              qualityThreshold=linkerIdentity[[x]], 
                              returnUnmatched=TRUE, returnLowScored=TRUE, 
                              doRC=doRC, ...) 
          }
        }        
      }
      names(linkerTrimmed) <- samplesToProcess
      
      ## check if any error occured during alignments ##
      if(any(grepl("simpleError",linkerTrimmed))) {
        stop("Error encountered in Linker Trimming functions",
             paste(names(linkerTrimmed[grepl("simpleError",linkerTrimmed)]),
                   collapse=", "))
      }
      
      ## remove samples with no linker hits from further processing ##
      culprits <- grep("No hits found",linkerTrimmed)
      if(length(culprits)>0) {
        message("Following sample(s) had no hits for Linker alignment: ",
                paste(samplesToProcess[culprits],collapse=", "))
        samplesToProcess <- samplesToProcess[-c(culprits)]
        linkerTrimmed <- linkerTrimmed[-c(culprits)]
      }        
      
      cleanit <- gc()
      
      if(!bypassChecks | showStats) {
        eval(expression(trimmedObj <- "linkerTrimmed", rawObj <- "toProcess", 
                        featureTrim <- "Linker", valueColname <- "PercOfDecoded"))
        .showFindStats()
      }
      
      message("Adding linker info back to the object")
      ## modify metadata attribute and write back to sampleInfo object
      ## for primerID based samples...write back all the returned attibutes
      for(x in names(linkerTrimmed)) {
        cat(".") 
        if(isPaired[[x]]) {
          for(y in unique(unlist(sapply(linkerTrimmed[[x]], names)))) {
            bore <- sapply(linkerTrimmed[[x]], "[[", y)
            newAttrName <- paste0(ifelse(y=="hits","",y),"linkered")
            sampleInfo$sectors[[sector]]$samples[[x]][[newAttrName]] <- bore
            rm(bore)
          }
        } else {
          for(y in names(linkerTrimmed[[x]])) {
            newAttrName <- paste0(ifelse(y=="hits","",y),"linkered")
            sampleInfo$sectors[[sector]]$samples[[x]][[newAttrName]] <- 
              linkerTrimmed[[x]][[y]]
          }
        }
      }      
      rm(linkerTrimmed)
      cleanit <- gc()
    }
  }
  sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
  return(sampleInfo)
}

#' Compare LTRed/Primed sequences to all linkers. 
#'
#' Given a SampleInfo object, the function compares LTRed sequences from each sample per sector to all the linker sequences present in the run. The output is a summary table of counts of good matches to all the linkers per sample. 
#'
#' @param sampleInfo sample information SimpleList object outputted from \code{\link{findPrimers}} or \code{\link{findLTRs}}, which holds decoded sequences for samples per sector/quadrant along with information of sample to primer associations.
#' @param qualityThreshold percent of linker length to match, round(nchar(linker)*qualityThreshold). Default is 0.55. Only applies to non-primerID based linkers
#' @param qualityThreshold1 percent of first part of patternSeq to match. Default is 0.75. Only applies to primerID based linker search.
#' @param qualityThreshold2 percent of second part of patternSeq to match. Default is 0.50. Only applies to primerID based linker search.
#' @param doRC perform reverse complement search of the linker sequence. Default is TRUE. Highly recommended!
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param samplenames a vector of samplenames to process. Default is NULL, which processes all samples from sampleInfo object.
#' @param ... extra parameters to be passed to \code{\link{pairwiseAlignment}}.
#'
#' @return a dataframe of counts. 
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}}, \code{\link{findPrimers}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #troubleshootLinkers(sampleInfo,showStats=TRUE)
#'
troubleshootLinkers <- function(sampleInfo, qualityThreshold=0.55, 
                                qualityThreshold1=0.75, qualityThreshold2=0.50, 
                                doRC=TRUE, parallel=TRUE, samplenames=NULL, ...) {    
  
  .checkArgs_SEQed()
  
  ## test if there are decoded sequences in the sampleinfo object ##
  toProcess <- extractFeature(sampleInfo, feature="decoded")
  toProcessSamples <- sapply(toProcess, names, simplify=FALSE)
  sectors <- names(toProcess)   
  rm(toProcess)
  cleanit <- gc()
  
  results <- data.frame()
  for(sector in sectors) {
    message("Processing sector ",sector)
    
    ## prepare sample to linker association ##
    sampleLinkers <- toupper(extractFeature(sampleInfo, sector=sector,
                                            feature="linkersequence")[[sector]])
    if(length(sampleLinkers)==0 | mean(nchar(sampleLinkers))<=10 | 
         all(is.na(sampleLinkers))) {
      stop("Either Linker sequence is too short (<=10) or ",
           "no Linkers found in sample information object.")
    }
    
    samplesToProcess <- toProcessSamples[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }
    
    toProcess.seqs <- extractSeqs(sampleInfo,sector,samplesToProcess,
                                  feature="decoded")[[sector]]    
    
    ## find paired end samples
    isPaired <- extractFeature(sampleInfo, sector, feature='pairedend')[[sector]]
      
    ## do all by all comparison of linkers ##
    message("\tFinding Linkers.")
    for(linkerSeq in unique(as.character(sampleLinkers))) {
      message("Checking ",linkerSeq)
      linkerTrimmed <- foreach(x=iter(samplesToProcess), .inorder=TRUE, 
                               .errorhandling="pass", 
                               .export=c("linkerSeq", "toProcess.seqs", "doRC",
                                         "qualityThreshold", "qualityThreshold1",
                                         "qualityThreshold2", "pairwiseAlignSeqs",
                                         "primerIDAlignSeqs", "isPaired"), 
                               .packages="Biostrings") %dopar% {
        if(length(unlist(gregexpr("N",linkerSeq)))>3) {
          if(isPaired[[x]]) {
            p1 <- try_default(length(primerIDAlignSeqs(toProcess.seqs[[x]]$pair1, 
                                                       linkerSeq, 
                                                       qualityThreshold1, 
                                                       qualityThreshold2,
                                                       doRC=doRC, ...)$hits), 
                              0, quiet=TRUE)

            p2 <- try_default(length(primerIDAlignSeqs(toProcess.seqs[[x]]$pair2, 
                                                       linkerSeq, 
                                                       qualityThreshold1, 
                                                       qualityThreshold2,
                                                       doRC=doRC, ...)$hits), 
                              0, quiet=TRUE)
            
            list("pair1"=p1, "pair2"=p2)
          } else {
            try_default(length(primerIDAlignSeqs(toProcess.seqs[[x]], linkerSeq, 
                                                 qualityThreshold1, qualityThreshold2, 
                                                 doRC=doRC, ...)$hits), 0, quiet=TRUE)
          }
        } else {
          if(isPaired[[x]]) {
            p1 <- try_default(length(pairwiseAlignSeqs(toProcess.seqs[[x]]$pair1,
                                                       linkerSeq, "middle",
                                                       qualityThreshold, doRC=doRC,
                                                       ...)), 0, quiet=TRUE)

            p2 <- try_default(length(pairwiseAlignSeqs(toProcess.seqs[[x]]$pair2,
                                                       linkerSeq, "right",
                                                       qualityThreshold, doRC=doRC,
                                                       ...)), 0, quiet=TRUE)
            
            list("pair1"=p1, "pair2"=p2)                       
          } else {
            try_default(length(pairwiseAlignSeqs(toProcess.seqs[[x]], linkerSeq,
                                                 "middle", qualityThreshold, 
                                                 doRC=doRC, ...)), 0, quiet=TRUE)
          }
        }
      }
      names(linkerTrimmed) <- samplesToProcess
      
      isPaired <- sapply(linkerTrimmed, is.list)
      if(any(isPaired)) {
        totalSeqs <- sapply(toProcess.seqs[isPaired], sapply, length)
        linkerhits <- sapply(linkerTrimmed[isPaired], unlist)
        PercentOfTotal <- linkerhits/totalSeqs[,names(linkerTrimmed[isPaired])]        
      } else {
        totalSeqs <- sapply(toProcess.seqs[!isPaired], length)
        linkerhits <- as.numeric(unlist(linkerTrimmed[!isPaired]))
        PercentOfTotal<- linkerhits/totalSeqs[names(linkerTrimmed[!isPaired])]
      }
      
      bore <- data.frame("linkerSeq"=linkerSeq, "samplename"=names(linkerTrimmed), 
                         "linkerhits"=linkerhits, "PercentOfTotal"=PercentOfTotal,
                         stringsAsFactors=FALSE)
      rownames(bore) <- NULL
      results <- rbind(results, bore)
      
      cleanit <- gc()
    }        
  }  
  sampleLinkers <- extractFeature(sampleInfo, feature="linkersequence")
  names(sampleLinkers) <- NULL
  sampleLinkers <- unlist(sampleLinkers)
  linkersample <- as.data.frame(sampleLinkers)
  linkersample <- tapply(rownames(linkersample), linkersample$sampleLinkers, 
                         paste, collapse=",")
  
  results$CorrectLinker <- with(results,
                                sampleLinkers[as.character(samplename)]==as.character(linkerSeq))
  results$CorrectSample <- with(results, linkersample[as.character(linkerSeq)])
  return(results)
}

#' Find and trim a short pattern sequence from the subject. 
#'
#' This function facilitates finding and trimming of a short pattern sequence from a collection of subject sequences. The trimming is dictated by side parameter. For more information on the trimming process see the 'side' parameter documentation in \code{\link{trimSeqs}}. For information regarding the pattern alignment see the documentation for \code{\link{pairwiseAlignSeqs}}. This function is meant for aligning a short pattern onto large collection of subjects. If you are looking to align a vector sequence to subjects, then please use BLAT.
#'
#' @param patternSeq DNAString object or a sequence containing the query sequence to search.
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern.
#' @param side which side of the sequence to perform the search & trimming: left, right or middle. Default is 'left'.
#' @param offBy integer value dictating if the trimming base should be offset by X number of bases. Default is 0.
#' @param alignWay method to utilize for detecting the primers. One of following: "slow" (Default), "fast", or "blat". Fast, calls \code{\link{vpairwiseAlignSeqs}} and uses \code{\link{vpatternMatch}} at its core, which is less accurate with indels and mismatches but much faster. Slow, calls \code{\link{pairwiseAlignSeqs}} and uses \code{\link{pairwiseAlignment}} at its core, which is accurate with indels and mismatches but slower. Blat will use \code{\link{blatSeqs}}.
#' @param ... parameters to be passed to \code{\link{pairwiseAlignment}}, \code{\link{vpairwiseAlignSeqs}} or \code{\link{blatSeqs}} depending on which method is defined in 'alignWay' parameter.
#'
#' @return DNAStringSet object with pattern sequence removed from the subject sequences. 
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{extractFeature}}, \code{\link{extractSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findPrimers}}, \code{\link{findLinkers}}
#'
#' @export
#'
#' @examples 
#' #findAndTrimSeq(patternSeq="AGACCCTTTT",subjectSeqs=DNAStringSet(c("AGACCCTTTTGAGCAGCAT","AGACCCTTGGTCGACTCA","AGACCCTTTTGACGAGCTAG")), qualityThreshold=.85, doRC=F, side="left", offBy=1, alignWay = "slow")
#'
findAndTrimSeq <- function(patternSeq, subjectSeqs, side = "left", offBy = 0, 
                           alignWay = "slow", ...) {
  
  .checkArgs_SEQed()
  
  coords <- switch(alignWay,
                   fast = vpairwiseAlignSeqs(subjectSeqs, patternSeq, side, ...),
                   slow = pairwiseAlignSeqs(subjectSeqs, patternSeq, side, ...),
                   blat = blatSeqs(subjectSeqs, patternSeq, ...)
  )
  
  res <- trimSeqs(subjectSeqs, coords, side, offBy)
  if(removeSubjectNamesAfter) {
    names(res) <- NULL
  }
  res
}

#' Find and trim vector sequence from reads. 
#'
#' This function facilitates finding and trimming of long/short fragments of vector present in LM-PCR products. The algorithm looks for vector sequence present anywhere within the read and trims according longest contiguous match on either end of the read.
#'
#' @param reads DNAStringSet object containing sequences to be trimmed for vector.
#' @param Vector DNAString object containing vector sequence to be searched in reads.
#' @param minLength integer value dictating minimum length of trimmed sequences to return. Default is 10.
#' @param returnCoords return the coordinates of vector start-stop for the matching reads. Defaults to FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @return DNAStringSet object with Vector sequence removed from the reads. If returnCoords=TRUE, then a list of two named elements "hits" & "reads". The first element, "hits" is a GRanges object with properties of matched region and whether it was considered valid denoted by 'good.row'. The second element, "reads" is a DNAStringSet object with Vector sequence removed from the reads.
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{pslToRangedObject}}, \code{\link{blatSeqs}}, \code{\link{read.blast8}}, \code{\link{findAndTrimSeq}}
#' @export
#'
#' @examples 
#' #findAndRemoveVector(reads, Vector)
#'
findAndRemoveVector <- function(reads, Vector, minLength=10, 
                                returnCoords=FALSE, parallel=TRUE) {
  
  .checkArgs_SEQed()
  
  hits <- read.blast8(blatSeqs(query=reads, subject=Vector, parallel=parallel,
                               blatParameters=c(stepSize = 3, tileSize = 8,
                                                minIdentity=70, minScore=5,
                                                repMatch = 112312, out = "blast8")))
  hits$evalue <- NULL
  hits <- reduce(pslToRangedObject(hits, useTargetAsRef=FALSE, isblast8=T),
                 min.gapwidth=10)
  hits$qSize <- width(reads[as.character(seqnames(hits))])
  
  ## only consider hits that start or end with vector...others would be wierdos! ##
  tocheck <- start(hits) <= 10 | hits$qSize-end(hits) <= 11
  hits <- hits[tocheck]
  
  ## queries where genomic sequence is surrounded by vector sequence (counts>1)...
  ## we will fix these later!
  qNames <- as.character(seqnames(hits))
  counts <- table(qNames)
  good.row <- !qNames %in% names(counts[counts>1])
  widthSum <- tapply(width(hits), qNames, sum)
  rm(counts)
  
  ### queries composed of recombined forms of vector
  qSizes <- hits$qSize
  allCovered <- qSizes==widthSum[qNames] | qSizes==widthSum[qNames]-1 | 
    qSizes==widthSum[qNames]+1 | qSizes==widthSum[qNames]+2 | 
    qSizes==widthSum[qNames]-2 | qSizes==widthSum[qNames]+3 | 
    qSizes==widthSum[qNames]-3
  rm("qSizes","widthSum")
  
  ### seqs w/ vector sequences with allCovered removed
  toTrim <- reads[names(reads) %in% qNames[!allCovered]]
  
  ### seqs w/o any vector sequences
  reads <- reads[!names(reads) %in% qNames] 
  
  ### queries where genomic sequence is either at the beginning or the end
  hits$qEndDiff <- hits$qSize-end(hits)
  hits$StartDiff <- start(hits)-0
  hits2 <- hits[good.row]
  
  # take the side which has the higher unmatched bases to the vector sequence #
  hits2$trueStart <- ifelse(hits2$StartDiff > hits2$qEndDiff, 0, end(hits2)+1)
  hits2$trueEnd <- ifelse(hits2$StartDiff > hits2$qEndDiff, 
                          start(hits2)-1, hits2$qSize)
  hits2$good.row <- !hits2$StartDiff==hits2$qEndDiff
  hits2$type <- "genomic either side"
  
  ### fix queries where genomic sequence is surrounded by vector sequence
  if(any(!good.row)) {
    hits.list <- split(hits[!good.row], qNames[!good.row])
    bores <- foreach(bore = iter(hits.list), .inorder = FALSE, 
                     .packages = c("GenomicRanges")) %dopar% {
                       trueRange <- gaps(ranges(bore))
                       bore$trueStart <- start(trueRange[width(trueRange)==
                                                           max(width(trueRange))])
                       bore$trueEnd <- end(trueRange[width(trueRange)==
                                                       max(width(trueRange))])
                       bore$type <- "genomic in middle"
                       bore$good.row <- TRUE
                       bore[1,names(mcols(hits2))]
                     }
    hits2 <- c(hits2, do.call(c,bores))
    rm("hits.list","bores")
  }
  rm("hits")
  
  atStart <- hits2$StartDiff<=10
  atEnd <- hits2$qEndDiff<=10
  
  ### remove sequences where vector sequence is in the middle 
  hits2$good.row[!atStart & !atEnd] <- FALSE 
  hits2$type[!atStart & !atEnd] <- "vector in middle" 
  good.row <- hits2$good.row
  
  ## extract sequence with trueStart & trueEnd 
  trimmed <- subseq(toTrim[as.character(seqnames(hits2))[good.row]],
                    start=hits2$trueStart[good.row]+1,
                    end=hits2$trueEnd[good.row])
  
  ### combined seqs.good with trimmed ###
  reads <- c(reads, trimmed[width(trimmed)>=minLength])
  
  if(returnCoords) {
    list("hits"=hits2, "reads"=reads)
  } else {
    reads
  }
}

#' Trim sequences from a specific side.
#'
#' This function trims a DNAStringSet object using the ranges from left, right, or middle of the sequence. This is a helper function utilized in \code{\link{primerIDAlignSeqs}} and \code{\link{extractSeqs}}. If dnaSet and coords are not the same length, then they are required to have a names attribute to perform the matched trimming. 
#'
#' @param dnaSet DNAStringSet object containing sequences to be trimmed.
#' @param coords IRanges object containing boundaries.
#' @param side either 'left','right',or the Default 'middle'.
#' @param offBy integer value dictating if the supplied coordinates should be offset by X number of bases. Default is 0.
#'
#' @return a DNAStringSet object with trimmed sequences. 
#'
#' @note If side is left, then any sequence following end of coords+offBy is returned. If side is right, then sequence preceding start of coords-offBy is returned. If side is middle, then sequence contained in coords is returned where offBy is added to start and subtracted from end in coords.
#'
#' @seealso \code{\link{extractSeqs}}, \code{\link{primerIDAlignSeqs}}
#'
#' @export
#'
#' @examples 
#' #trimSeqs(dnaSet,coords)
#' #trimSeqs(dnaSet,coords,side="left",offBy=1)
#'
trimSeqs <- function(dnaSet, coords, side="middle", offBy=0) {
  stopifnot(class(dnaSet) %in% c("DNAStringSet", "DNAString"))
  stopifnot(class(coords)=="IRanges")
  
  if(length(dnaSet)==0 | length(coords)==0) {
    stop("dnaSet/coords is empty. Please supply reads/coords to be trimmed.")
  }
  
  # check if both dnaSet and coords has 'names' attribute, 
  # if yes then check if they have matching names, else check lengths. 
  if(is.null(names(dnaSet)) | is.null(names(coords))) {
    stopifnot(length(dnaSet)==length(coords))
  } else {
    rows <- match(names(coords), names(dnaSet))
    if(any(is.na(rows))) {
      stop("Some of the names in coords are not present in dnaSet")
    }
    if(!is.ordered(rows)) {
      dnaSet <- dnaSet[rows]
      if(!identical(names(dnaSet), names(coords))) {
        stop("Names are not identical between dnaSet and coords parameter")
      }
    }
  }        
  
  # temp helper function to show messages #
  .showMessage <- function(x) {
    message("Following sequences were removed from trimming since their ",
            "coordinates+offBy were out of sequence length: ", 
            paste(x,collapse=", "))
  }
  
  # trim by side and check if any of the coords are off the sequence length in dnaSet
  if(tolower(side)=="left") {
    test <- end(coords)+offBy > width(dnaSet) | end(coords)+offBy < 1
    if(any(test)) {
      .showMessage(names(dnaSet)[test])
    }
    subseq(dnaSet[!test], start=end(coords[!test])+offBy)
  } else if (tolower(side)=="right") {
    test <- start(coords)-offBy > width(dnaSet) | end(coords)+offBy < 1
    if(any(test)) {
      .showMessage(names(dnaSet)[test])
    }
    subseq(dnaSet[!test], end=start(coords[!test])-offBy)
  } else {
    test <- start(coords)+offBy > width(dnaSet) | 
      end(coords)-offBy > width(dnaSet) | start(coords)+offBy < 1
    if(any(test)) {
      .showMessage(names(dnaSet)[test])
    }
    subseq(dnaSet[!test], 
           start=start(coords[!test])+offBy, 
           end=end(coords[!test])-offBy)
  }    
}

#' Extract sequences for a feature in the sampleInfo object.
#'
#' Given a sampleInfo object, the function extracts sequences for a defined feature.
#'
#' @param sampleInfo sample information SimpleList object, which samples per sector/quadrant information along with other metadata.
#' @param sector specific sector to extract sequences from. Default is NULL, which extracts all sectors. 
#' @param samplename specific sample to extract sequences from. Default is NULL, which extracts all samples. 
#' @param feature which part of sequence to extract (case sensitive). Options include: primed, !primed, LTRed, !LTRed, linkered, !linkered, primerIDs, genomic, genomicLinkered, decoded, and unDecoded. If a sample was primerIDed and processed by \code{\link{primerIDAlignSeqs}}, then all the rejected and unmatched attributes can be prepended to the feature. Example: Rejectedlinkered, RejectedprimerIDslinkered, Absentlinkered, or unAnchoredprimerIDslinkered. When feature is genomic, it includes sequences which are primed, LTRed, linkered, and !linkered. The genomicLinkered is same as genomic minus the !linkered. When feature is decoded, it includes everything that demultiplexed. The '!' in front of a feature extracts the inverse. One can only get unDecoded sequences if returnUnmatched was TRUE in \code{\link{decodeByBarcode}}.
#' @param trim whether to trim the given feature from sequences or keep it. Default is TRUE. This option is ignored for feature with '!'.
#' @param minReadLength threshold for minimum length of trimmed sequences to return.
#' @param sideReturn if trim=TRUE, which side of the sequence to return: left, middle, or right. Defaults to NULL and determined automatically. Doesn't apply to features: decoded, genomic or genomicLinkered.
#' @param pairReturn if the data is paired end, then from which pair to return the feature. Options are "pair1", "pair2", or defaults to "both". Ignored if data is single end. 
#'
#' @return a listed DNAStringSet object structed by sector then sample. Note: when feature='genomic' or 'genomicLinkered' and when data is paired end, then "pair2" includes union of reads from both pairs which found LTR.
#'
#' @seealso \code{\link{findPrimers}}, \code{\link{findLTRs}}, \code{\link{findLinkers}}, \code{\link{trimSeqs}}, \code{\link{extractFeature}}, \code{\link{getSectorsForSamples}}
#'
#' @export
#'
#' @examples 
#' #extractSeqs(sampleInfo)
#' #extractSeqs(sampleInfo,feature="primed")
#'
extractSeqs <- function(sampleInfo, sector=NULL, samplename=NULL, feature="genomic",
                        trim=TRUE, minReadLength=1, sideReturn=NULL, 
                        pairReturn="both") {
  
  .checkArgs_SEQed()
  
  # get all sectors and samplenames in each sector or a specific sector
  res <- getSectorsForSamples(sampleInfo, sector, samplename)
  sectors <- res[["sectors"]]
  samplenames <- res[["samplenames"]]
  
  if(feature=="unDecoded") {
    sapply(sectors, function(y) { 
      allmetadata <- metadata(sampleInfo$sectors[[y]])
      if("unDecodedSeqs" %in% names(allmetadata)) {
        allmetadata$unDecodedSeqs 
      } else {
        message("No unDecoded attribute found for the supplied ",
                "sampleInfo object and sector ", y)
      }
    })
  } else {        
    res <- sapply(sectors, function(y) {
      sapply(samplenames[[y]], function(x,y) { 
        decoded <- sampleInfo$sectors[[y]]$samples[[x]]$decoded
        isPaired <- sampleInfo$sectors[[y]]$samples[[x]]$pairedend
        if(isPaired & pairReturn!='both') {
          decoded <- decoded[[pairReturn]]
        }

        if(feature!="decoded") {
          # get ride of ! from feature, else R wont know what to do when making a new object with ! in the front.
          assign(gsub("!","",feature), 
                 sampleInfo$sectors[[y]]$samples[[x]][[gsub("!","",feature)]])
        }
        
        if(feature=="decoded") {
          decoded
        } else if (feature %in% c("genomic","genomicLinkered")) {
          primed <- sampleInfo$sectors[[y]]$samples[[x]]$primed
          LTRed <- sampleInfo$sectors[[y]]$samples[[x]]$LTRed
          linkered <- sampleInfo$sectors[[y]]$samples[[x]]$linkered
          
          if(isPaired & pairReturn!='both') {
            primed <- primed[[pairReturn]]
            LTRed <- LTRed[[pairReturn]]
            linkered <- linkered[[pairReturn]]
          }
          
          # if reads dont have LTRs...use primers instead...but throw a warning!
          LTRed.test <- if(is.list(LTRed)) {
            any(sapply(LTRed,is.null))
          } else {
            is.null(LTRed)
          }
          if(LTRed.test) {
            warning("LTRed information not found for",x,
                    "...using primer end as starting boundary.", immediate.=TRUE)                                
          }
          
          # test if there are any primed and/or linkered reads
          p.l.test <- if(is.list(primed)) {
            any(sapply(primed,is.null), sapply(linkered,is.null))
          } else {
            any(is.null(primed), is.null(linkered))
          }
          
          if(p.l.test) { 
            message("No sequences found for requested feature (", feature,
                    ") for sample: ", x,"...skipping.") 
          } else {
            if(trim) {
              # get all LTRed reads and make ends = size of each read
              if(isPaired & pairReturn=='both') {
                if(LTRed.test) {
                  LTRed <- primed
                }

                starts <- lapply(LTRed, function(p) structure(end(p), 
                                                              names=names(p)))
                ## add reads present pair1 to pair2 where LTR wasn't found
                ## make things consistent to pair1 atleast!
                loners <- setdiff(names(starts$pair1), names(starts$pair2))
                starts$pair2 <- c(starts$pair2, 
                                  structure(rep(0L,length(loners)), names=loners))
                ends <- lapply(decoded, function(p) structure(width(p), 
                                                              names=names(p)))

                stopifnot(identical(names(starts), names(ends)))
                stopifnot(mapply(function(s,e) all(names(s) %in% names(e)),
                                 starts, ends))

                coords <- mapply(function(s,e) {
                  IRanges(start=s+1, end=e[names(s)], names=names(s))
                }, starts, ends)

                # alter ends for reads where linker was present
                # fix cases where start of linker earlier than edge of LTR/primer end
                ends <- lapply(linkered, function(p) structure(start(p), 
                                                               names=names(p)))
                stopifnot(identical(names(coords), names(ends)))
                
                coords <- mapply(function(p, e) {
                  there <- intersect(names(p), names(e))
                  e <- as.numeric(e[there]-1)
                  culprits <- start(p[there]) > e
                  end(p[there][!culprits]) <- e[!culprits]
                  end(p[there][culprits]) <- start(p[there][culprits])
                  p
                }, coords, ends)

                if(feature=="genomicLinkered") { 
                  stopifnot(identical(names(coords), names(linkered)))
                  coords <- mapply(function(p, l) p[names(p) %in% names(l)], 
                                   coords, linkered)
                }
                
                # trim it and return non zero length sequences
                if(any(sapply(coords, length)>0)) {
                  stopifnot(identical(names(coords), names(decoded)))
                  mapply(function(d, p) {
                    seqs <- trimSeqs(d,p)
                    seqs[width(seqs)>=minReadLength]
                  }, decoded, coords)
                } else {
                  message("No linkered reads found for sample: ",x,"...skipping.")
                }
              } else {
                if(LTRed.test) {
                  LTRed <- primed
                }
                
                starts <- structure(end(LTRed), names=names(LTRed))
                ends <- structure(width(decoded), names=names(decoded))

                stopifnot(identical(names(starts), names(ends[names(starts)])))
                coords <- IRanges(start=starts+1,
                                  end=ends[names(starts)],
                                  names=names(starts))
                
                # alter ends for reads where linker was present
                # fix cases where start of linker earlier than edge of LTR/primer end
                ends <- structure(start(linkered), names=names(linkered))
                there <- intersect(names(coords), names(ends))
                ends <- as.numeric(ends[there]-1)
                culprits <- start(coords[there]) > ends
                end(coords[there][!culprits]) <- ends[!culprits]
                end(coords[there][culprits]) <- start(coords[there][culprits])
                                
                if(feature=="genomicLinkered") { 
                  coords <- coords[names(coords) %in% names(linkered)] 
                }
                
                # trim it and return non zero length sequences
                if(length(coords)>0) {
                  seqs <- trimSeqs(decoded,coords)
                  seqs[width(seqs)>=minReadLength]
                } else {
                  message("No linkered reads found for sample: ",x,"...skipping.")
                }                
              }
            } else {
              # just retuning ranges...simply match names from decoded to the request
              if(isPaired & pairReturn=='both') {
                if(feature=="genomicLinkered") { 
                  mapply(function(d,l) d[names(d) %in% names(l)],
                         decoded, linkered)
                } else { ## everything past primer and LTR
                  if(LTRed.test) {
                    mapply(function(d,l) d[names(d) %in% names(l)],
                           decoded, primed)
                  } else {
                    mapply(function(d,l) d[names(d) %in% names(l)],
                           decoded, LTRed)
                  }
                }
              } else {
                if(feature=="genomicLinkered") { 
                  decoded[names(decoded) %in% names(linkered)]
                } else { ## everything past primer and LTR
                  if(is.null(LTRed)) {
                    decoded[names(decoded) %in% names(primed)]
                  } else {
                    decoded[names(decoded) %in% names(LTRed)]
                  }
                }
              }              
            }
          }
        } else {
          notFeature <- grepl("!",feature)
          ## no need to trim if looking for "not" based feature since there are no coordinates for it
          if(notFeature) {
            if(isPaired & pairReturn=='both') {
              stopifnot(identical(names(decoded), names(get(gsub("!","",feature)))))
              toreturn <- mapply(function(d,g) d[!names(d) %in% names(g)], 
                                 decoded, get(gsub("!","",feature)))
              
              toreturn <- list()
              toreturn[["pair1"]] <- decoded$pair1[!names(decoded$pair1) %in% 
                                                     names(get(gsub("!","",
                                                                    feature))$pair1)]
              toreturn[["pair2"]] <- decoded$pair2[!names(decoded$pair2) %in% 
                                                     names(get(gsub("!","",
                                                                    feature))$pair2)]
            } else {
              toreturn <- decoded[!names(decoded) %in% 
                                    names(get(gsub("!","",feature)))]
            }
 
            if(length(toreturn)) {
              toreturn
            } else {
              message("No sequences found for requested feature (", feature,
                      ") for sample: ", x,"...skipping.") 
            }
          } else {                    
            if(is.null(get(feature))) { 
              message("No sequences found for requested feature (", feature,
                      ") for sample: ", x,"...skipping.") 
            } else {
              if(isPaired & pairReturn=='both') {
                stopifnot(identical(names(decoded), names(get(feature))))                
                res.seq <- list("pair1"=decoded$pair1[names(decoded$pair1) %in% 
                                                        names(get(feature)$pair1)],
                                "pair2"=decoded$pair2[names(decoded$pair2) %in% 
                                                        names(get(feature)$pair2)])
              } else {
                res.seq <- decoded[names(decoded) %in% names(get(feature))]  
              }
              
              if(trim) {
                if(is.null(sideReturn)) {
                  sidetype <- ifelse(grepl("primerID",feature,ignore.case=TRUE), 
                                     "middle", 
                                     ifelse(grepl("primed|LTRed",feature,
                                                  ignore.case=TRUE), 
                                            "left", 
                                            ifelse(grepl("linkered",feature,
                                                         ignore.case=TRUE), 
                                                   "right", 
                                                   "middle")))
                } else {
                  sidetype <- tolower(sideReturn)
                }                            
                
                offByLength <- ifelse(sidetype=="middle",0,1)
                
                if(isPaired & pairReturn=='both') {
                  seqs.p1 <- trimSeqs(res.seq$pair1, get(feature)$pair1, 
                                      side=sidetype, offBy=offByLength)
                  seqs.p1 <- seqs.p1[width(seqs.p1)>=minReadLength]
                  
                  seqs.p2 <- trimSeqs(res.seq$pair2, get(feature)$pair2, 
                                      side=sidetype, offBy=offByLength)
                  seqs.p2 <- seqs.p2[width(seqs.p2)>=minReadLength]
                  
                  list("pair1"=seqs.p1, "pair2"=seqs.p2)
                } else {
                  seqs <- trimSeqs(res.seq, get(feature), 
                                   side=sidetype, offBy=offByLength)
                  seqs[width(seqs)>=minReadLength]
                }                
              } else {
                res.seq
              }
            }
          }
        }
      }, y=y)
    }, simplify=FALSE)
    
    simpletons <- !sapply(res,class)=="matrix"
    if(any(simpletons)) {
      lengthTest <- lapply(lapply(lapply(res[simpletons], 
                                         function(x) sapply(x, length)),">",0),
                           which)
      res <- mapply(function(x,y) x[y], res[simpletons], lengthTest, SIMPLIFY=FALSE)  
    }
    
    if(any(!simpletons)) {
      res[!simpletons] <- sapply(res[!simpletons], as, 'DataFrame')
    }
    
    res
  }
}

#' Extract a specific feature/attribute of the sampleInfo object.
#'
#' Given a sampleInfo object, the function extracts a defined feature(s) for given sample or sector.
#'
#' @param sampleInfo sample information SimpleList object, which samples per sector/quadrant information along with other metadata.
#' @param sector a vector or specific sector to extract the feature from. Default is NULL, which extracts all sectors. 
#' @param samplename a character vector or a specific sample to extract feature from. Default is NULL, which extracts all samples. 
#' @param feature Options include: primed, LTRed, linkered, decoded, and any of the metadata. Default is NULL. When feature='metadata', then it returns names of all the metadata elements associated with the sample as a comma separated list.
#'
#' @return a list or list of lists depending upon which parameters were supplied.
#'
#' @seealso \code{\link{addFeature}}, \code{\link{findPrimers}}, \code{\link{findLTRs}}, \code{\link{findLinkers}}, \code{\link{extractSeqs}}, \code{\link{trimSeqs}}, \code{\link{getSectorsForSamples}}
#'
#' @export
#'
#' @examples 
#' #extractFeature(sampleInfo,feature="primed")
#'
extractFeature <- function(sampleInfo, sector=NULL, samplename=NULL, feature=NULL) {
  
  .checkArgs_SEQed()
  
  # get all sectors and samplenames in each sector or a specific sector
  res <- getSectorsForSamples(sampleInfo, sector, samplename)
  sectors <- res[["sectors"]]
  samplenames <- res[["samplenames"]]
  
  res <- sapply(sectors, function(y) {
    sapply(samplenames[[y]], function(x,y) { 
      if(feature=="metadata") {
        paste(names(sampleInfo$sectors[[y]]$samples[[x]]), collapse=", ")
      } else {
        res <- sampleInfo$sectors[[y]]$samples[[x]][[feature]]
        ## convert any factor based vector to the appropriate regular vector
        if(class(res)=="factor") { 
          if(!any(is.na(suppressWarnings(as.numeric(levels(res)))))) { 
            as.numeric(as.character(res)) } 
          else { 
            as.character(res) 
          }    
        } else if (class(res)=="character"){ 
          ## if a numeric vector is stored as character, convert it back to numeric
          if(!any(is.na(suppressWarnings(as.numeric(res))))) { 
            as.numeric(res) 
          } else { 
            res 
          }
        } else {
          res
        }
      }
    }, y=y)
  }, simplify=FALSE)
  
  simpletons <- !sapply(res,class)=="matrix"
  if(any(simpletons)) {
    lengthTest <- lapply(lapply(lapply(res[simpletons], 
                                       function(x) sapply(x, length)),">",0),
                         which)
    res <- mapply(function(x,y) x[y], res[simpletons], lengthTest, SIMPLIFY=FALSE)  
  }
  
  if(any(!simpletons)) {
    res[!simpletons] <- sapply(res[!simpletons], as, 'DataFrame')
  }
  
  res
}

#' Add a specific feature/attribute to the sampleInfo object.
#'
#' Given a sampleInfo object, the function adds a new feature for the given samples & sectors.
#'
#' @param sampleInfo sample information SimpleList object, which samples per sector/quadrant information along with other metadata.
#' @param sector a vector or a specific sector to add the new feature(s) to. Default is NULL, in which case the sectors are searched from samplename parameter.
#' @param samplename a character vector or a specific sample to add the new feature(s) to. Default is NULL.
#' @param feature a string of naming the new feature to add for the defined samplename and sector.
#' @param value the value or a named list of samplenames & values which is assigned for the defined sector, samplename, and feature. Example: list("Sample1"="ACDTDASD")
#'
#' @return modified sampleInfo object with new feature(s) added.
#'
#' @seealso \code{\link{findPrimers}}, \code{\link{extractSeqs}}, \code{\link{trimSeqs}}, \code{\link{extractFeature}}, \code{\link{getSectorsForSamples}}
#'
#' @export
#'
#' @examples 
#' #addFeature(sampleInfo,"1","Sample1",feature="primerIDed","ACDTDASD")
#'
addFeature <- function(sampleInfo, sector=NULL, samplename=NULL, feature=NULL, 
                       value=NULL) {
  
  .checkArgs_SEQed()
  
  if(is.null(value)) {
    stop("Please define the value parameter.")
  }
  
  if(is.null(samplename)) {
    stop("Please define the samplename to add the features to.")
  }
  
  if(!all(samplename %in% names(value))) {
    stop("Not all samplename(s) are found in value parameter")
  }
  
  # get all sectors and samplenames in each sector or a specific sector
  res <- getSectorsForSamples(sampleInfo, sector, samplename)
  sectors <- res[["sectors"]]
  samplenames <- res[["samplenames"]]
  
  for(y in sectors) {
    for(x in samplenames[[y]]) {
      cat(".")
      sampleInfo$sectors[[y]]$samples[[x]][[feature]] <- value[[x]]
    }
  }
  return(sampleInfo)
}

#' Get sectors for samples defined in the sampleInfo object.
#'
#' Given a sampleInfo object, the function gets the sectors for each samplename. This is an accessory function utilized by other functions of this package to aid sector retrieval.
#'
#' @param sampleInfo sample information SimpleList object, which samples per sector/quadrant information along with other metadata.
#' @param sector a specific sector or vector of sectors if known ahead of time. Default is NULL, which extracts all sectors. 
#' @param samplename a specific sample or vector of samplenames to get sectors for. Default is NULL, which extracts all samples. 
#' @param returnDf return results in a dataframe. Default is FALSE.
#'
#' @return If returnDf=TRUE, then a dataframe of sector associated with each samplename, else a named list of length two: x[["sectors"]] and x[["samplenames"]]
#'
#' @seealso \code{\link{extractSeqs}}, \code{\link{extractFeature}}, \code{\link{addFeature}}
#'
#' @export
#'
#' @examples 
#' #getSectorsForSamples(sampleInfo,samplename="SampleName1")
#'
getSectorsForSamples <- function(sampleInfo, sector=NULL, samplename=NULL,
                                 returnDf=FALSE) {
  
  .checkArgs_SEQed()
  
  # get all sectors and samplenames in each sector or a specific sector if defined
  if(is.null(sector)) {
    sectors <- names(sampleInfo$sectors)        
  } else {
    sectors <- sector
  }
  samplenames <- sapply(sectors,
                        function(x) names(sampleInfo$sectors[[x]]$samples),
                        simplify=FALSE)
  
  # if specific samplename(s) is defined, then search to find where they are
  if(is.null(samplename)) {
    samplename <- unlist(samplenames)     
  }
  
  if(!all(samplename %in% unlist(samplenames))) {            
    stop("Following sample(s) do not exist on given sector(s) (",
         paste(sectors,collapse=", "),") in the supplied sampleInfo object: ",
         paste(samplename[!samplename %in% unlist(samplenames)],collapse=", "))
  }
  sectors <- names(which(unlist(lapply(lapply(samplenames,"%in%",samplename), any)))) 
  if(returnDf) {
    return(do.call(rbind, lapply(sectors, function(x) { 
      data.frame(samplename=samplenames[[x]][samplenames[[x]] %in% samplename], 
                 sector=x, stringsAsFactors=FALSE)
    })))
  } else {
    samplenames <- sapply(sectors, 
                          function(x) 
                            samplenames[[x]][samplenames[[x]] %in% samplename],
                          simplify=FALSE)
    return(list("sectors"=sectors, "samplenames"=samplenames))
  }
}

#' Read fasta/fastq given the path or sampleInfo object.
#'
#' Given a sequence reads file path, the function returns a DNAStringSet object.
#'
#' @param seqFilePath a path to fasta/fastq reads file or a sampleInfo object returned by \code{\link{read.SeqFolder}}
#' @param sector specific sector to reads sequences from. Default is 1, and not required if seqFilePath is a direct file path rather than sampleInfo object.
#'  @param isPaired does the sector contain paired end reads? Default is FALSE
#'
#' @return if isPaired is FALSE, then a DNAStringSet object, else a list of DNAStringSet objects of three elements corresponding to reads from "barcode", "pair1", and "pair2". Note: "pair2" is reverse complemented!
#'
#' @seealso \code{\link{decodeByBarcode}}, \code{\link{read.SeqFolder}}, \code{\link{extractSeqs}}
#'
#' @export
#'
read.seqsFromSector <- function(seqFilePath=NULL, sector=1, isPaired=FALSE) {
  if(is.null(seqFilePath)) {
    stop("Missing seqFilePath!")
  }
  
  if(class(seqFilePath)=="SimpleList") {
    seqfilePattern <- seqFilePath$seqfilePattern
    seqFilePaths <- seqFilePath$seqFilePaths
    
    if(isPaired) {
      filePath <- normalizePath(grep(paste0(gsub("R1|R2|I1",".*",sector),
                                            seqfilePattern), seqFilePaths,
                                     value=TRUE), mustWork=TRUE)
      
      if(length(filePath)==0) {
        stop("No sequence file found for sector: ", sector,
             " in seqFilePath variable (",
             paste(seqFilePaths,collapse=" \n"), ") using pattern (",
             paste0(sector,seqfilePattern),")")
      }
      
      if(!any(grepl("I1", basename(filePath)))) {
        warning("No index/barcode file (I1) found.")
      }
      
      pair2 <- paste0(sub("(.*)R\\d.*","\\1",sector),
                      ifelse(sub(".*(R\\d).*","\\1",sector)=="R2","R1","R2"))
      if(!any(grepl(pair2, basename(filePath)))) {
        stop("Pair #2 or Linker end read file not present!")
      }
      seqFilePath <- filePath
    
    } else {
      filePath <- normalizePath(grep(paste0(sector,seqfilePattern),
                                     seqFilePaths, value=TRUE), mustWork=TRUE)
      
      if(length(filePath)==0) {
        stop("No sequence file found for sector: ", sector,
             " in seqFilePath variable (",
             paste(seqFilePaths,collapse=" \n"), ") using pattern (",
             paste0(sector,seqfilePattern),")")
      }
      
      if(length(filePath)>1) {
        stop("Multiple sequence file found for sector: ", sector,
             " in seqFilePath variable (", paste(seqFilePaths,collapse=" \n"),
             ") using pattern (", paste0(sector,seqfilePattern),")")
      }
      seqFilePath <- filePath
    }
  }
  
  message("Reading:\n", paste(seqFilePath, collapse="\n"))
  if(any(grepl("fastq",seqFilePath,ignore.case=TRUE))) {    
    dnaSet <- sapply(seqFilePath, function(x) {
      bore <- readFastq(x)
      dnaSet <- sread(bore)
      names(dnaSet) <- sub("^\\S+-(\\S+) .+$", "\\1", id(bore), perl=TRUE)
      if(any(duplicated(names(dnaSet)))) {
        stop("Duplicate definition lines found in file: ", x)
      }
      if(length(dnaSet)==0) {
        stop("No sequences found in file: ", x)
      }
      dnaSet
    })
        
    if(isPaired & length(seqFilePath)>1) {
      LTRside <- grep(sector, basename(names(dnaSet)))
      #names(dnaSet[[LTRside]]) <- paste0("@pair1side@",names(dnaSet[[LTRside]]))
      
      linkerSide <- grep(pair2, basename(names(dnaSet)))      
      #names(dnaSet[[linkerSide]]) <- paste0("@pair2side@",names(dnaSet[[linkerSide]]))
      
      barcodes <- grep("I1", basename(names(dnaSet)))
      
      dnaSet <- list("barcode"=dnaSet[[barcodes]],
                     "pair1"=dnaSet[[LTRside]],
                     "pair2"=reverseComplement(dnaSet[[linkerSide]]))
    } else {
      ## for single-end data!
      dnaSet <- dnaSet[[1]]
    }
  } else {
    dnaSet <- readDNAStringSet(seqFilePath)
    if(any(duplicated(names(dnaSet)))) {
      stop("Duplicate definition lines found in file(s): ", 
           paste(seqFilePath, collapse=" * "))
    }
    
    if(length(dnaSet)==0) {
      stop("No sequences found")
    }
  }
  
  return(dnaSet)
}

#' Write a fasta file per sample in parallel
#'
#' Given a listed DNAStringSet object return from \code{\link{extractSeqs}}, the function writes a fasta file for each sample as defined in filePath parameter.
#'
#' @param dnaSet listed DNAStringSet object containing sequences to be written.
#' @param filePath a path write the fasta files per sample. Default is current working directory.
#' @param filePrefix prefix the filenames with a string. Default is 'processed' followed by samplename.
#' @param prependSamplenames Prepend definition lines with samplenames. Default is TRUE. Make sure the dnaSet parameter is a named list where names are used as samplenames.
#' @param format either fasta (the default) or fastq.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @seealso \code{\link{decodeByBarcode}}, \code{\link{read.SeqFolder}}, \code{\link{extractSeqs}}
#'
#' @note
#' \itemize{
#'   \item Writing of the files is done using \code{\link{writeXStringSet}} with parameter append=TRUE. This is to aggregate reads from a sample which might be present in more than one sector. 
#'   \item If data is paired end, then each pair will be written separately with designations in the filename.
#'   \item If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#' }
#'
#' @export
#'
#' @examples 
#' #write.listedDNAStringSet(dnaSet)
#'
write.listedDNAStringSet <- function(dnaSet, filePath=".", filePrefix="processed", 
                                     prependSamplenames=TRUE, format="fasta", 
                                     parallel=FALSE) {
  stopifnot(class(dnaSet)=="list")
  
  if(!parallel) { registerDoSEQ() }
  if(filePrefix=="") { filePrefix <- NA }
  
  for(x in seq_along(dnaSet)) {
    foreach(outputSeqs=iter(as(dnaSet[[x]],"list")), 
            samplename=iter(names(dnaSet[[x]])),
            .inorder=TRUE, .errorhandling="stop", .packages="Biostrings", 
            .export=c("filePath", "filePrefix", "prependSamplenames")) %dopar% {
      
      if(is.list(outputSeqs)) {
        list.test <- any(sapply(outputSeqs, length)>0)
      } else {
        list.test <- length(outputSeqs)>0
        outputSeqs <- list("tempy"=outputSeqs)
      }
      
      if(list.test) {
        for(p in names(outputSeqs)) {
          pairname <- ifelse(p=="tempy", NA, p)
          
          if(is.null(names(outputSeqs[[p]]))) {
            message("No names attribute found for ", samplename,
                    " ... using artifically generated names")
            names(outputSeqs[[p]]) <- paste("read",1:length(outputSeqs[[p]]),sep="-")
          }
          
          ## remove '.' at the beginning of the filename incase filePrefix is empty
          filename <- paste(na.omit(c(filePrefix, samplename, pairname, 
                                      ifelse(format=="fastq", "fastq", "fa"))), 
                            collapse=".")
          filename <- paste(filePath, filename, sep="/")
          
          if(prependSamplenames) {
            names(outputSeqs[[p]]) <- paste(samplename, names(outputSeqs[[p]]), sep="-")
          }
          
          writeXStringSet(outputSeqs[[p]], file=filename, format=format, append=TRUE) 
        } 
      } else {
        message("No reads written for ", samplename)
      } 
    }
  }
}

#' Check args and set defaults for functions dealing with reads.
#'
#' This function checks all the arguments passed to a function related to aligning or trimming a read and then sets default values for internal use. Evaluation of this function happens in the parent function.
#'
.checkArgs_SEQed <- function() {
  
  checks <- expression( 
    if("subjectSeqs" %in% names(formals())) { 
      if(is.null(subjectSeqs) | length(subjectSeqs)==0) {
        stop("subjectSeqs paramter is empty. Please supply reads to be aligned")
      }      
      
      ## give names if not there for troubleshooting purpose in later steps
      removeSubjectNamesAfter <- FALSE
      if(is.null(names(subjectSeqs))) {
        removeSubjectNamesAfter <- TRUE
        names(subjectSeqs) <- paste("read", 1:length(subjectSeqs))
      } 
    },
    
    if("patternSeq" %in% names(formals())) { 
      if(is.null(patternSeq) | length(patternSeq)==0) {
        stop("patternSeq paramter is empty. Please supply reads to be aligned")
      } else if (length(patternSeq)>1) {
        stop("More than 1 patternSeq is defined. Please only supply one pattern.")
      }
    },
    
    if("reads" %in% names(formals())) { 
      if(is.null(reads) | length(reads)==0) {
        stop("reads paramter is empty. Please supply reads to be aligned")
      }
    },
    
    if("Vector" %in% names(formals())) { 
      if(is.null(Vector) | length(Vector)==0) {
        stop("Vector paramter is empty. Please supply Vector to be aligned")
      }
    },
    
    if("parallel" %in% names(formals())) { 
      if(!parallel) { 
        registerDoSEQ() 
      }
    },
    
    if("sampleInfo" %in% names(formals())) { 
      stopifnot(is(sampleInfo,"SimpleList"))
    },

    if("feature" %in% names(formals())) {
      if(is.null(feature)) {
        stop("Please define a feature to extract.")
      }
    }      
  )
  
  eval.parent(checks)
}

#' Check args and set defaults for functions dealing with hits of aligned reads.
#'
#' This function checks all the arguments passed to a function related to checking hits of an aligned read and then sets default values for internal use. Evaluation of this function happens in the parent function.
#'
.checkArgsSetDefaults_ALIGNed <- function() {
  
  checks <- expression()
  
  eval.parent(checks)
}

#' Start a gfServer instance
#'
#' Start a gfServer indexed reference genome to align batch of sequences using BLAT gfServer/gfClient protocol.
#'
#' @param seqDir absolute or relative path to the genome index (nib/2bit files).
#' @param host name of the machine to run gfServer on. Default: localhost
#' @param port a port number to host the gfServer with. Default is 5560.
#' @param gfServerOpts a character vector of options to be passed to gfServer command on top of server defaults. Default: c(repMatch=112312, stepSize=5, tileSize=10). Set this to NULL to start gfServer with defaults.
#'
#' @return system command status for executing gfServer command.
#'
#' @seealso \code{\link{stopgfServer}}, \code{\link{read.psl}}, \code{\link{blatSeqs}}, \code{\link{read.blast8}}
#'
#' @export
#'
#' @examples 
#' #startgfServer(seqDir="/usr/local/blatSuite34/hg18.2bit",port=5560)
#' #stopgfServer(port=5560)
#' 
startgfServer <- function(seqDir=NULL, host="localhost", port=5560, 
                          gfServerOpts=c(repMatch=112312, stepSize=5, tileSize=10)) {
  if(is.null(seqDir)) {
    stop("Please define the path of nib/2bit files containing the ",
         "indexed reference sequence(s)")
  }
  
  cmd <- sprintf("gfServer start %s %i %s %s &", 
                 host, port, 
                 ifelse(!is.null(gfServerOpts),
                        paste(paste("-",names(gfServerOpts),sep=""), 
                              gfServerOpts, collapse=" ", sep="="),
                        ""), 
                 normalizePath(seqDir))
  message(cmd)
  system(cmd)        
  
  ## wait for server to load & be ready
  message("Loading BLAT server...please wait.")    
  searchCMD <- sprintf("gfServer status %s %s", host, port)
  while(system(searchCMD,ignore.stderr=TRUE)!=0) {
    cat(".")
    Sys.sleep(10)
  }
}

#' Stop a gfServer instance
#'
#' Take down the gfServer instance.
#'
#' @param host name of the machine running gfServer. Default: localhost
#' @param port the same number you started the gfServer with. 
#'
#' @return system command status for executing kill command.
#'
#' @seealso \code{\link{startgfServer}}
#'
#' @export
#'
#' @examples 
#' 	#stopgfServer(port=5560)
#' 
stopgfServer <- function(host="localhost", port=NULL) {
  if(is.null(port)) {
    stop("Please define the port gfServer is running on.")
  }
  
  cmd <- sprintf("kill `ps ax | grep '%s' | grep -v 'grep' | awk '{print $1}'`", 
                 paste("gfServer start", host, port))
  system(cmd)
}

#' Align a listed DNAStringSet to a reference using gfClient or standalone BLAT.
#'
#' Align sequences from a listed DNAStringSet object returned from \code{\link{extractSeqs}} to an indexed reference genome using gfServer/gfClient protocol or using standalone BLAT and return the psl file as a GRanges object. This function heavily relies on defaults of \code{\link{blatSeqs}}.
#'
#' @param dnaSetList DNAStringSet object containing sequences to be aligned against the reference.
#' @param ... parameters to be passed to \code{\link{blatSeqs}}.
#'
#' @return a list of GRanges object reflecting psl file type per set of sequences.
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{startgfServer}}, \code{\link{stopgfServer}}, \code{\link{blatSeqs}}, \code{\link{read.psl}}, \code{\link{pslToRangedObject}}, \code{\link{read.blast8}}
#'
#' @export
#' 
blatListedSet <- function(dnaSetList=NULL, ...) {
  if(is.null(dnaSetList)) {
    stop("dnaSetList is empty. Please supply a listed DNAStringSet ",
         "returned from extractSeqs() to be aligned against a reference")
  }
  
  sapply(names(dnaSetList), function(x) {
    do.call(GRangesList,sapply(names(dnaSetList[[x]]), function(y) {
      if(length(dnaSetList[[x]][[y]])>0) {
        message("BLATing ",y)
        outFiles <- blatSeqs(query=dnaSetList[[x]][[y]], ...)
        read.psl(outFiles, bestScoring=TRUE, 
                 asGRanges=TRUE, removeFile=TRUE, parallel=FALSE)
      }
    }))
  })
}

#' Convert psl dataframe to RangedData/GRanges
#'
#' Convert psl dataframe to RangedData or GRanges object using either the query or target as the reference data column. 
#'
#' @param x dataframe reflecting psl format
#' @param useTargetAsRef use target(tName) or query(qName) as the chromosome or the reference data. Default is TRUE.
#' @param asGRanges make a GRanges object instead of RangedData. Default is TRUE.
#' @param isblast8 the input dataframe blast8 format output from BLAT. Default is FALSE.
#'
#' @return a GRanges object reflecting psl file type. If asGRanges=FALSE, then RangedData object is returned.
#'
#' @seealso \code{\link{read.psl}}, \code{\link{read.blast8}}, \code{\link{blatListedSet}}
#'
#' @export
#'
#' @examples 
#' #pslToRangedObject(psl)
#' #pslToRangedObject(psl, asGRanges=FALSE)
#' #pslToRangedObject(psl, useTargetAsRef=FALSE)
#'
pslToRangedObject <- function(x, useTargetAsRef=TRUE, asGRanges=TRUE, isblast8=FALSE) {
  if(useTargetAsRef) {
    metadataCols <- c(setdiff(names(x), c("tName","tStart","tEnd","strand")),
                      ifelse(isblast8, NA, "tStarts"))
    out <- GRanges(seqnames=x$tName, IRanges(start=x$tStart, end=x$tEnd),
                   strand=x$strand)     
  } else {
    metadataCols <- c(setdiff(names(x), c("qName","qStart","qEnd","strand")),
                      ifelse(isblast8, NA, "qStarts"))
    out <- GRanges(seqnames=x$qName, IRanges(start=x$qStart, end=x$qEnd),
                   strand=x$strand)
  }
  
  for(f in na.omit(metadataCols)) {
    mcols(out)[[f]] <- x[,f]
  } 
  
  if(!asGRanges) {
    out <- as(out, "RangedData")
  }
  
  out
}

#' Split DNA sequences into smaller files.
#'
#' Given a vector of sequences or DNAStringSet or a FASTA filename, the function splits it into smaller pieces as denoted by totalFiles parameter.
#'
#' @param x a DNAStringSet object, or a FASTA filename.
#' @param totalFiles an integer indicating how many files to create. Default is 4.
#' @param suffix a word to add to each file created. Default is "tempy".
#' @param filename name of the file if x is a DNAStringSet object. Default is "queryFile.fa".
#'
#' @return a vector of filename names created.
#'
#' @seealso \code{\link{blatSeqs}}
#'
#' @export
#'
#' @examples 
#' #splitSeqsToFiles(dnaSeq,10,"tempyQ","myDNAseqs.fa")
#' #splitSeqsToFiles("mySeqs.fa",5,"tempyQ")
#'
splitSeqsToFiles <- function(x, totalFiles=4, suffix="tempy", 
                             filename="queryFile.fa") {
  if(is.atomic(x)) {
    message("Splitting file ",x)
    totalSeqs <- length(fasta.info(x, use.names=FALSE))
    chunks <- round(totalSeqs/totalFiles)
    ## incase totalSeqs is lower than number of files to be created!
    chunks <- ifelse(chunks>0, chunks, totalSeqs) 
    
    starts <- seq(0,totalSeqs,by=chunks) ## create chunks of starts
    for(skippy in starts[starts!=totalSeqs]) {
      filename.out <- paste(x, skippy, suffix,sep=".")
      ## no need to read the entire file...save memory by reading in N lines
      query.tmp <- readBStringSet(x,nrec=chunks, skip=skippy) 
      writeXStringSet(query.tmp, file=filename.out, format="fasta")            
    }
    return(list.files(path=dirname(x), 
                      pattern=paste0(basename(x),".*", suffix, "$"), 
                      full.names=TRUE))
  } else if (class(x)=="DNAStringSet") {
    message("Splitting Reads.")
    totalSeqs <- length(x)
    chunks <- round(totalSeqs/totalFiles)
    starts <- seq(1, totalSeqs, by=chunks)
    stops <- unique(c(seq(chunks, totalSeqs, by=chunks), totalSeqs))
    stopifnot(length(starts)==length(stops))        
    for(skippy in 1:length(starts)) {
      filename.out <- paste(filename, skippy, suffix,sep=".")            
      writeXStringSet(x[starts[skippy]:stops[skippy]], file=filename.out,
                      format="fasta")            
    }            
    return(list.files(path=".", 
                      pattern=paste(filename,".*",suffix,"$",sep=""), 
                      full.names=TRUE))
  } else {
    stop("Dont know what is supplied in parameter x.")
  }
}

#' Align sequences using BLAT.
#'
#' Align batch of sequences using standalone BLAT or gfServer/gfClient protocol for alignment against an indexed reference genome. Depending on parameters provided, the function either aligns batch of files to a reference genome using gfClient or takes sequences from query & subject parameters and aligns them using standalone BLAT. If standaloneBlat=FALSE and gfServer is not launched apriori, this function will start one using \code{\link{startgfServer}} and kill it using \code{\link{stopgfServer}} upon successful execution. 
#'
#' @param query an object of DNAStringSet, a character vector of filename(s), or a path/pattern of fasta files to BLAT. Default is NULL.
#' @param subject an object of DNAStringSet, a character vector, or a path to an indexed genome (nibs,2bits) to serve as a reference or target to the query. Default is NULL. If the subject is a path to a nib or 2bit file, then standaloneBlat will not work!
#' @param standaloneBlat use standalone BLAT as suppose to gfServer/gfClient protocol. Default is TRUE.
#' @param port the same number you started the gfServer with. Required if standaloneBlat=FALSE. Default is 5560.
#' @param host name of the machine running gfServer. Default is 'localhost' and only used when standaloneBlat=FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param gzipResults gzip the output files? Default is TRUE.
#' @param blatParameters a character vector of options to be passed to gfClient/BLAT command except for 'nohead' option. Default: c(minIdentity=70, minScore=5, stepSize=5, tileSize=10, repMatch=112312, dots=50, q="dna", t="dna", out="psl"). Be sure to only pass parameters accepted by either BLAT or gfClient. For example, if repMatch or stepSize parameters are specified when using gfClient, then the function will simply ignore them! The defaults are configured to align a 19bp sequence with 70\% identity.
#'
#' @return a character vector of psl filenames.
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{startgfServer}}, \code{\link{stopgfServer}}, \code{\link{read.psl}}, \code{\link{splitSeqsToFiles}}, \code{\link{read.blast8}}
#'
#' @export
#'
#' @examples 
#' #blatSeqs(dnaSeqs, subjectSeqs, blatParameters=c(minIdentity=90, minScore=10, tileSize=10, dots=10, q="dna", t="dna", out="blast8"))
#' #blatSeqs(dnaSeqs, "/usr/local/genomeIndex/hg18.2bit", standaloneBlat=FALSE)
#' #blatSeqs("mySeqs.fa", "/usr/local/genomeIndex/hg18.2bit", standaloneBlat=FALSE)
#' #' #blatSeqs("my.*.fa", "/usr/local/genomeIndex/hg18.2bit", standaloneBlat=FALSE)
#'
blatSeqs <- function(query=NULL, subject=NULL, standaloneBlat=TRUE, port=5560, 
                     host="localhost", parallel=TRUE, gzipResults=TRUE, 
                     blatParameters=c(minIdentity=70, minScore=5, stepSize=5, 
                                      tileSize=10, repMatch=112312, dots=50, 
                                      q="dna", t="dna", out="psl")) {
  
  ## get all BLAT options from the system for comparison to blatParameters later
  suppressWarnings(blatOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", 
                                          grep("\\s+-.+=", 
                                               system("blat",intern=TRUE), 
                                               value=TRUE))))
  
  suppressWarnings(gfClientOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", 
                                              grep("\\s+-.+=", 
                                                   system("gfClient", intern=TRUE), 
                                                   value=TRUE))))
                                                   
  gfServerOpts <- c("tileSize","stepSize","minMatch","maxGap", "trans","log",
                    "seqLog","syslog","logFacility","mask","repMatch","maxDnaHits",
                    "maxTransHits","maxNtSize","maxAsSize","canStop")                                               
  
  if(!standaloneBlat) { 
    message("Using gfClient protocol to perform BLAT.")
    if(is.null(port)) {
      stop("The port paramter is empty. ",
           "Please define the port used to start gfServer with")
    }
  }
  
  ## check the subject parameter
  if(is.null(subject) | length(subject)==0) {
    stop("The subject parameter is empty. ", 
         "Please supply subject sequences or a path to 2bit or nib files to ",
         "serve as reference/target")
  } else {
    subjectFile <- NULL
    if(is.atomic(subject)) {
      if (any(grepl("\\.2bit$|\\.nib$", subject))) {
        if(standaloneBlat) { 
          stop("Standalone BLAT cannot be used when subject is an indexed ",
               "nib or 2bit file.") 
        }
        indexFileDir <- dirname(subject)
        subjectFile <- list.files(path=indexFileDir, 
                                  pattern=basename(subject), full.names=TRUE)
        if(length(subjectFile)==0) { 
          stop("The file(s) supplied in subject parameter doesn't exist.") }
      } else {
        ## change object type if necessary for troubleshooting purpose in later steps
        subject <- DNAStringSet(subject)
      }
    }
    
    if(is.null(subjectFile)) {
      ## subjectFile is still null so it means that subject is a DNAStringSet
      if(is.null(names(subject))) { ## add names of subject if not present
        names(subject) <- paste("subject", 1:length(subject))
      }
      
      ## write out the subject sequences into a fasta file
      filename.seq <- "subjectFile.fa.tempyS"
      writeXStringSet(subject, file=filename.seq, format="fasta")                                  
      subjectFile <- filename.seq
    }
  }
  
  ## check the query parameter
  if(is.null(query) | length(query)==0) {
    stop("The query parameter is empty. Please supply reads to be aligned")
  } else {
    queryFiles <- NULL
    if(is.atomic(query)) {
      if (any(grepl("\\.fna$|\\.fa$|\\*", query))) {
        ## detect whether query paramter is a regex or list of files
        if(any(grepl("\\*|\\$|\\+|\\^",query))) {
          queryFiles <- list.files(path=dirname(query), pattern=basename(query), 
                                   full.names=TRUE)            
        } else {
          queryFiles <- query
        }
        
        if(parallel) {
          ## split the fasta files into smaller chunks for parallel BLATing
          queryFiles <- unlist(sapply(queryFiles,
                                      function(f) splitSeqsToFiles(f, 
                                                                   getDoParWorkers(),
                                                                   "tempyQ")), 
                               use.names=FALSE)                    
        }
      } else {
        ## change object type if necessary for troubleshooting purpose in later steps
        query <- DNAStringSet(query)
      }
    }
    
    if(is.null(queryFiles)) {
      ## queryFiles is still null so it means that query is a DNAStringSet           
      if(is.null(names(query))) {  ## fix names of query if not present
        names(query) <- paste("read", 1:length(query),sep="-")
      }  
      
      ## write out the query sequences into fasta files
      if(parallel) {
        queryFiles <- splitSeqsToFiles(query, getDoParWorkers(), "tempyQ")
      } else {
        queryFiles <- "queryFile.fa.tempyQ"
        writeXStringSet(query, file=queryFiles, format="fasta")                
      }
    }
  }
  
  ## perform the Blatting of queryFiles vs subjectFile/indexFiles using gfClient/standalone BLAT
  if(!parallel) { registerDoSEQ() }
  
  ## do some formatting ##
  queryFiles <- as.character(queryFiles)
  subjectFile <- as.character(subjectFile)
  
  ## BLAT it ##
  if(standaloneBlat) {        
    blatOpts <- blatParameters[names(blatParameters) %in% blatOpts]
    stopifnot(length(subjectFile)==1)
    filenames <- foreach(x=iter(queryFiles), .inorder=FALSE,
                         .export=c("blatOpts","subjectFile","gzipResults")) %dopar% {
      filename.out <- paste(x, blatOpts["out"], sep=".")
      cmd <- paste("blat", paste(paste("-",names(blatOpts),sep=""), blatOpts, 
                                 collapse=" ", sep="="), "-noHead", 
                   subjectFile, x, filename.out)
      message(cmd)
      system(cmd)
      
      ## no need to save splitted files!
      if(grepl("\\.tempyQ$",x)) { system(paste("rm",x)) } 
      
      if(gzipResults) { 
        system(paste("gzip", filename.out))
        filename.out <- paste(filename.out, "gz", sep=".") 
      }
      filename.out
    }
    if(grepl("\\.tempyS$",subjectFile)) { system(paste("rm", subjectFile)) }
  } else {
    # start the gfServer if not started already! #
    killFlag <- FALSE
    searchCMD <- sprintf("gfServer status %s %s", host, port)
    if(system(searchCMD,ignore.stderr=TRUE)!=0) {
      message("Starting gfServer.")        
      startgfServer(seqDir=subjectFile, host=host, port=port, 
                    gfServerOpts=blatParameters[names(blatParameters) %in% gfServerOpts])
      killFlag <- TRUE
    }         
    
    gfClientOpts <- blatParameters[names(blatParameters) %in% gfClientOpts]
    stopifnot(length(subjectFile)>0)
    filenames <- foreach(x=iter(queryFiles), .inorder=FALSE, 
                         .export=c("gfClientOpts","host","port",
                                   "indexFileDir","gzipResults")) %dopar% {
      filename.out <- paste(x, gfClientOpts["out"], sep=".")
      cmd <- paste("gfClient",
                   paste(paste("-",names(gfClientOpts),sep=""), gfClientOpts, 
                         collapse=" ", sep="="), "-nohead", 
                   host, port, "/", x, filename.out)
      message(cmd)
      system(cmd)
      
      ## no need to save splitted files!
      if(grepl("\\.tempyQ$", x)) { system(paste("rm", x)) } 
      
      if(gzipResults) { 
        system(paste("gzip", filename.out))
        filename.out <- paste(filename.out, "gz", sep=".") 
      }
      filename.out
    }  
    
    ## only kill if the gfServer was started from within this function ## 
    if(killFlag) {
      # stop to conserve memory #
      message("Kill gfServer.")        
      stopgfServer(port=port)      
    }
  }
  return(unlist(filenames))
}

#' Read psl file(s) outputted by BLAT
#'
#' Given filename(s), the function reads the psl file format from BLAT as a data frame and performs basic score filtering if indicated. Any other file format will yield errors or erroneous results. Make sure there is no header row!
#'
#' @param pslFile psl filename, or vector of filenames, or a pattern of files to import.
#' @param bestScoring report only best scoring hits instead of all hits. Default is TRUE. Score is calculated by matches-misMatches-qBaseInsert-tBaseInsert.
#' @param asRangedData return a RangedData object instead of a dataframe. Default is FALSE
#' @param asGRanges return a GRanges object instead of a dataframe. Default is FALSE
#' @param removeFile remove the psl file(s) after importing. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @return a dataframe reflecting psl file type. If asGRanges=T or asRangedData=T, then a GRanges object or RangedData object, respectively.  
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{startgfServer}}, \code{\link{blatSeqs}}, \code{\link{read.blast8}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #read.psl(pslFile="processed.*.psl$")
#' #read.psl(pslFile=c("sample1hits.psl","sample2hits.psl"))
#'
read.psl <- function(pslFile=NULL, bestScoring=TRUE, asRangedData=FALSE, 
                     asGRanges=FALSE, removeFile=TRUE, parallel=FALSE) {
  if(is.null(pslFile) | length(pslFile)==0) {
    stop("pslFile parameter empty. Please supply a filename to be read.")
  }
  
  if (any(grepl("\\*",pslFile))) {
    ## vector of filenames
    pslFile <- list.files(path=dirname(pslFile), 
                          pattern=basename(pslFile), full.names=TRUE)      
  }
  
  if(length(pslFile)==0) { 
    stop("No file(s) found with given paramter in pslFile:", pslFile) 
  }
  
  if(!parallel) { registerDoSEQ() }
  
  ## setup psl columns + classes
  cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", 
            "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", 
            "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", 
            "blockSizes", "qStarts", "tStarts")
  cols.class <- c(rep("numeric",8), rep("character",2), rep("numeric",3),
                  "character", rep("numeric",4), rep("character",3))
  
  hits <- foreach(x=iter(pslFile), .inorder=FALSE, 
                  .export=c("cols","cols.class",
                            "bestScoring","removeFile")) %dopar% {
    message(x)
    hits.temp <- read.delim(x, header=FALSE, col.names=cols, 
                            stringsAsFactors=FALSE, colClasses=cols.class)
    if(removeFile) { system(paste("rm", x)) }
    if(bestScoring) {  
      ## do round one of bestScore here to reduce file size          
      hits.temp$score <- with(hits.temp, matches-misMatches-qBaseInsert-tBaseInsert)
      hits.temp <- arrange(hits.temp, qName, desc(score))
      bestScore <- with(hits.temp, tapply(score, qName, max))
      isBest <- with(hits.temp, score==bestScore[qName])
      hits.temp <- hits.temp[isBest,]
      rm("isBest","bestScore")
    }
    hits.temp
  }
  hits <- unique(do.call(rbind, hits))
  
  if(nrow(hits)==0) {
    stop("No hits found")
  }
  
  message("Ordering by qName")
  hits <- arrange(hits, qName)
  
  ## do round two of bestScore incase any got missed in round one
  if(bestScoring) {
    message("\t cherry picking!")
    hits$score <- with(hits, matches-misMatches-qBaseInsert-tBaseInsert)
    hits <- arrange(hits, qName, desc(score))
    bestScore <- with(hits, tapply(score,qName,max))
    isBest <- with(hits, score==bestScore[qName])
    hits <- hits[isBest,]
    rm("isBest","bestScore")
  }
  
  if(asRangedData | asGRanges) {
    hits <- pslToRangedObject(hits, useTargetAsRef=TRUE, asGRanges=asGRanges)
  }
  
  return(hits)
}

#' Read blast8 file(s) outputted by BLAT
#'
#' Given filename(s), the function reads the blast8 file format from BLAT as a data frame and performs basic score filtering if indicated. Any other file format will yield errors or erroneous results.
#'
#' @param files blast8 filename, or vector of filenames, or a pattern of files to import.
#' @param asRangedData return a RangedData object instead of a dataframe. Default is TRUE. Saves memory!
#' @param asGRanges return a GRanges object instead of a dataframe. Default is FALSE
#' @param removeFile remove the blast8 file(s) after importing. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @return a dataframe or RangedData object reflecting blast8 file type.
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{startgfServer}}, \code{\link{blatSeqs}}
#'
#' @export
#'
#' @examples 
#' #read.blast8(files="processed.*.blast8$")
#' #read.blast8(files=c("sample1hits.blast8","sample2hits.blast8"))
#'
read.blast8 <- function(files=NULL, asRangedData=FALSE, asGRanges=FALSE,
                        removeFile=TRUE, parallel=FALSE) {
  if(is.null(files) | length(files)==0) {
    stop("files parameter empty. Please supply a filename to be read.")
  }
  
  if (any(grepl("\\*",files))) {
    ## vector of filenames
    files <- list.files(path=dirname(files), 
                        pattern=basename(files), full.names=TRUE)
  }
  
  if(length(files)==0) { 
    stop("No file(s) found with given paramter in files:", files) 
  }
  
  if(!parallel) { registerDoSEQ() }
  
  ## setup blast8 columns + classes
  cols <- c("qName", "tName", "identity", "span", "misMatches", "gaps", "qStart", 
            "qEnd", "tStart", "tEnd", "evalue", "bitscore")
  cols.class<-c(rep("character",2), rep("numeric",10))
  
  hits <- foreach(x=iter(files), .inorder=FALSE, 
                  .export=c("cols","cols.class",
                            "bestScoring","removeFile")) %dopar% {
    message(x)
    hits.temp <- read.delim(x, header=FALSE, col.names=cols, 
                            stringsAsFactors=FALSE, colClasses=cols.class)
    hits.temp$strand <- with(hits.temp, ifelse(tStart>tEnd,"-","+"))
    # switch tStart & tEnd for cases where strand=='-' since it's reversed in blast8 format.
    rows <- hits.temp$strand=='-'
    tstarts <- hits.temp$tEnd[rows]
    tends <- hits.temp$tStart[rows]
    hits.temp$tStart[rows] <- tstarts
    hits.temp$tEnd[rows] <- tends
    rm("tstarts","tends","rows")
    
    if(removeFile) { system(paste("rm", x)) }
    hits.temp
  }
  hits <- unique(do.call(rbind, hits))
  
  if(nrow(hits)==0) {
    stop("No hits found")
  }
  
  message("Ordering by qName")
  hits <- arrange(hits,qName)
  
  if(asRangedData | asGRanges) {
    hits <- pslToRangedObject(hits, useTargetAsRef=TRUE, 
                              isblast8=TRUE, asGRanges=asGRanges)
  }
  
  return(hits)
}

#' Obtain integration sites from BLAT output
#'
#' Given a RangedData/GRanges object from \code{\link{read.psl}}, the function uses specified filtering parameters to obtain integration sites and maintain sequence attrition. The function will remove any non-best scoring alignments from the object if not already filtered apriori.
#'
#' @param psl.rd a RangedData/GRanges object reflecting psl format where tName is the spaces/seqnames.
#' @param startWithin upper bound limit of where the alignment should start within the query. Default is 3.
#' @param alignRatioThreshold cuttoff for (alignment span/read length). Default is 0.7.
#' @param genomicPercentIdentity cuttoff for (1-(misMatches/matches)). Default is 0.98.
#' @param correctByqStart use qStart to correct genomic position. This would account for sequencing/trimming errors. Position=ifelse(strand=="+",tStart-qStart,tEnd+qStart). Default is TRUE.
#' @param oneBased the coordinates in psl files are "zero based half open". The first base in a sequence is numbered zero rather than one. Enabling this would add +1 to the start and leave the end as is. Default is FALSE.
#'
#' @return a RangedData/GRanges object with integration sites which passed all filtering criteria. Each filtering parameter creates a new column to flag if a sequence/read passed that filter which follows the scheme: 'pass.FilterName'. Integration Site is marked by new column named 'Position'.
#'
#' @seealso \code{\link{startgfServer}}, \code{\link{read.psl}}, \code{\link{blatSeqs}}, \code{\link{blatListedSet}}, \code{\link{findIntegrations}}, \code{\link{pslToRangedObject}}, \code{\link{clusterSites}}, \code{\link{otuSites2}}, \code{\link{crossOverCheck}}, \code{\link{read.blast8}}
#'
#' @export
#'
#' @examples 
#' #getIntegrationSites(test.psl.rd)
#'
getIntegrationSites <- function(psl.rd=NULL, startWithin=3, alignRatioThreshold=0.7, 
                                genomicPercentIdentity=0.98, correctByqStart=TRUE, 
                                oneBased=FALSE) {
  stopifnot((is(psl.rd,"RangedData") | is(psl.rd,"GRanges")) & 
              !is.null(psl.rd) & !is.null(startWithin) & !length(psl.rd)==0 &
              !is.null(alignRatioThreshold) & !is.null(genomicPercentIdentity))
  
  isRangedData <- FALSE
  if(is(psl.rd, "RangedData")) {
    isRangedData <- TRUE
    psl.rd <- as(psl.rd, "GRanges")
  }
  
  ## check if required columns exist ##
  absentCols <- setdiff(c("qName", "qStart", "qSize","matches", "misMatches", 
                          "qBaseInsert", "tBaseInsert"), 
                        colnames(mcols(psl.rd)))
  						  
  if(length(absentCols)>0) {
  	stop("Following columns are absent from psl.rd object: ",
  		 paste(absentCols,collapse=","))
  }

  ## get the integration position by correcting for any insertions due to sequencing errors ##    
  if(correctByqStart) {
    mcols(psl.rd)$Position <- ifelse(as.character(strand(psl.rd))=="+",
                                     start(psl.rd)-mcols(psl.rd)$qStart,
                                     end(psl.rd)+mcols(psl.rd)$qStart)
  } else {
    mcols(psl.rd)$Position <- start(flank(psl.rd, width=-1))
  }
  
  ## get +1 based coordinate ##    
  if(oneBased) {
    mcols(psl.rd)$Position <- ifelse(as.character(strand(psl.rd))=="+",
                                     mcols(psl.rd)$Position+1,
                                     mcols(psl.rd)$Position)
  }
  
  # get scores for picking best hits and identify multihits later
  # check if scoring filtering hasn't already been applied by blat functions
  if(!"score" %in% colnames(mcols(psl.rd))) {
    message("Adding score column.")
    mcols(psl.rd)$score <- with(as.data.frame(mcols(psl.rd)), 
                                matches-misMatches-qBaseInsert-tBaseInsert)
    bestScore <- tapply(mcols(psl.rd)$score, as.character(mcols(psl.rd)$qName), max)
    isBest <- mcols(psl.rd)$score==bestScore[as.character(mcols(psl.rd)$qName)]
    psl.rd <- psl.rd[isBest,]
    rm("isBest","bestScore")
    cleanit <- gc()
  }    
  
  message("Performing QC checks.")
  # remove rows where the best hit dont start within first X bp
  mcols(psl.rd)$pass.startWithin <- mcols(psl.rd)$qStart<=startWithin
  
  # check if aligned ratio matches the threshold
  mcols(psl.rd)$alignRatio <- mcols(psl.rd)$score/mcols(psl.rd)$qSize
  mcols(psl.rd)$pass.alignRatio <- mcols(psl.rd)$alignRatio >= alignRatioThreshold
  
  # check for %identity    
  mcols(psl.rd)$percIdentity <- 1-(mcols(psl.rd)$misMatches/mcols(psl.rd)$matches)
  mcols(psl.rd)$pass.percIdentity <- mcols(psl.rd)$percIdentity >= genomicPercentIdentity
  
  ## find which query aligned to multiple places with equally good score aka...multihits
  cloneHits <- table(mcols(psl.rd)$qName)
  cloneHits <- cloneHits[as.character(mcols(psl.rd)$qName)]>1
  mcols(psl.rd)$isMultiHit <- as.logical(cloneHits)
  rm(cloneHits)    
  cleanit <- gc()
  
  mcols(psl.rd)$pass.allQC <- mcols(psl.rd)$pass.percIdentity & 
  								            mcols(psl.rd)$pass.alignRatio & 
  								            mcols(psl.rd)$pass.startWithin
  
  if(isRangedData) {
    psl.rd <- as(psl.rd, "RangedData")
  }
  
  return(psl.rd)
}

#' Cluster/Correct values within a window based on their frequency given discrete factors
#'
#' Given a group of discrete factors (i.e. position ids) and integer values, the function tries to correct/cluster the integer values based on their frequency in a defined windowsize.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand). Required if psl.rd parameter is not defined. 
#' @param value a vector of integer with values that needs to corrected/clustered (i.e. Positions). Required if psl.rd parameter is not defined. 
#' @param grouping additional vector of grouping by which to pool the rows (i.e. samplenames). Default is NULL.
#' @param psl.rd a RangedData/GRanges object returned from \code{\link{getIntegrationSites}}. Default is NULL. 
#' @param weight a numeric vector of weights to use when calculating frequency of value by posID and grouping if specified. Default is NULL.
#' @param windowSize size of window within which values should be corrected/clustered. Default is 5.
#' @param byQuartile flag denoting whether quartile based technique should be employed. See notes for details. Default is TRUE.
#' @param quartile if byQuartile=TRUE, then the quartile which serves as the threshold. Default is 0.70.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}. Process is split by the grouping the column.
#'
#' @note The algorithm for clustering when byQuartile=TRUE is as follows: for all values in each grouping, get a distribution and test if their frequency is >= quartile threshold. For values below the quartile threshold, test if any values overlap with the ones that passed the threshold and is within the defined windowSize. If there is a match, then merge with higher value, else leave it as is. This is only useful if the distribution is wide and polynodal. When byQuartile=FALSE, for each group the values within the defined window are merged with the next highest frequently occuring value, if freuquencies are tied then lowest value is used to represent the cluster. When psl.rd is passed, then multihits are ignored and only unique sites are clustered. All multihits will be tagged as a good 'clusterTopHit'.
#'
#' @return a data frame with clusteredValues and frequency shown alongside with the original input. If psl.rd parameter is defined then a RangedData/GRanges object is returned with three new columns appended at the end: clusteredPosition, clonecount, and clusterTopHit (a representative for a given cluster chosen by best scoring hit!). 
#'
#' @seealso \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, \code{\link{otuSites}}, \code{\link{otuSites2}}, \code{\link{crossOverCheck}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #clusterSites(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), value=c(rep(1000,2),5832,1000,12324,65738,928042), grouping=c('a','a','a','b','b','b','c'))
#' #clusterSites(grouping=test.psl.rd$grouping, psl.rd=test.psl.rd)
#'
clusterSites <- function(posID=NULL, value=NULL, grouping=NULL, psl.rd=NULL, 
                         weight=NULL, windowSize=5, byQuartile=FALSE, 
                         quartile=0.70, parallel=TRUE) {
  if(!parallel) { registerDoSEQ() }      
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(value))
  } else {
    
    isRangedData <- FALSE
    if(is(psl.rd, "RangedData")) {
      isRangedData <- TRUE
      psl.rd <- as(psl.rd, "GRanges")
    }
    
    ## dereplicate & converge sites! ##
    if(!"Position" %in% colnames(mcols(psl.rd))) {
      stop("The object supplied in psl.rd parameter does not have Position ",
           "column in it. Did you run getIntegrationSites() on it?")
    }
    
    ## make sure column(s) for picking best hits are there! ##
    isthere <- grepl("score", colnames(mcols(psl.rd)), ignore.case=TRUE)  	
    if(!any(isthere)) {
      message("No 'score' column found in the data. Using 'qEnd' as an ",
              "alternative to pick the best hit!")
      isthere <- grepl("qEnd", colnames(mcols(psl.rd)), ignore.case=TRUE)		
      if(!any(isthere)) {
        stop("No 'qEnd' column found in the data either...can't pick the ",
             "cluster best hit without 'qEnd' or 'score' column present in ",
             "the object supplied in psl.rd :(")
      }            
    }
    
    good.row <- rep(TRUE, length(psl.rd))
    multi.there <- grepl("isMultiHit", colnames(mcols(psl.rd)), ignore.case=TRUE)		
    if(any(multi.there)) { ## see if multihit column exists
      message("Found 'isMultiHit' column in the data. ",
              "These rows will be ignored for the calculation.")
      good.row <- good.row & !mcols(psl.rd)[[which(multi.there)]]
    }
    
    posIDs <- paste0(as.character(seqnames(psl.rd)), as.character(strand(psl.rd)))
    values <- mcols(psl.rd)$Position
    if(is.null(weight)) { ## see if sequences were dereplicated before in the pipeline which adds counts=x identifier to the deflines
      weight <- suppressWarnings(as.numeric(sub(".+counts=(\\d+)", "\\1",
                                                mcols(psl.rd)$qName)))
      if(all(is.na(weight))) { weight <- NULL } else { weight[is.na(weight)] <- 1 }
    }
    grouping <- if(is.null(grouping)) { rep("group1",length(values)) } else { grouping }
    clusters <- clusterSites(posIDs[good.row], 
                             values[good.row], 
                             grouping=grouping[good.row], 
                             weight=weight[good.row], 
                             windowSize=windowSize, 
                             byQuartile=byQuartile, quartile=quartile)
    
    message("Adding clustered data back to psl.rd.")        
    clusteredValues <- with(clusters,
                            split(clusteredValue, paste0(posID,value,grouping)))
    groupingVals <- paste0(posIDs, values, grouping)[good.row]
    mcols(psl.rd)$clusteredPosition <- mcols(psl.rd)$Position
    mcols(psl.rd)$clusteredPosition[good.row] <- as.numeric(clusteredValues[groupingVals])
    
    ## add frequency of new clusteredPosition ##
    clusteredValueFreq <- with(clusters, 
                               split(clusteredValue.freq, paste0(posID,value,grouping)))
    mcols(psl.rd)$clonecount <- 0
    mcols(psl.rd)$clonecount[good.row] <- as.numeric(clusteredValueFreq[groupingVals])
    rm("clusteredValueFreq","clusteredValues","clusters")
    cleanit <- gc()
    
    ## pick best scoring hit to represent a cluster ##
    ## make sure to avoid multihit rows! ##
    message("Picking best scoring hit to represent a cluster.")   
    groupingVals <- paste0(as.character(seqnames(psl.rd)), 
                           as.character(strand(psl.rd)), 
                           mcols(psl.rd)$clusteredPosition, grouping)
    bestScore <- tapply(mcols(psl.rd)[[which(isthere)]][good.row], 
                        groupingVals[good.row], max)
    isBest <- mcols(psl.rd)[[which(isthere)]][good.row] == 
      bestScore[groupingVals[good.row]]
    
    ## pick the first match for cases where >1 reads with the same coordinate had the same best scores ##
    tocheck <- which(isBest)
    res <- tapply(tocheck, names(tocheck), "[[", 1) 
    
    mcols(psl.rd)$clusterTopHit <- FALSE
    mcols(psl.rd)$clusterTopHit[good.row][res] <- TRUE
    mcols(psl.rd)$clusterTopHit[!good.row] <- TRUE
    
    message("Cleaning up!")
    rm("isBest","bestScore","posIDs","values","groupingVals")
    cleanit <- gc()
    
    if(isRangedData) {
      psl.rd <- as(psl.rd, "RangedData")
    }
    
    return(psl.rd)
  }
  
  # get frequencies of each posID & value combination by grouping #
  groups <- if(is.null(grouping)) { "" } else { grouping }
  weight2 <- if(is.null(weight)) { 1 } else { weight }
  sites <- arrange(data.frame(posID, value, grouping=groups, 
                              weight=weight2, posID2=paste0(groups, posID), 
                              stringsAsFactors=FALSE), posID2, value)
  sites <- count(sites, c("posID","value","grouping","posID2"), wt_var="weight")    
  rm("groups","weight2")
  
  if(byQuartile) {
    message("Clustering by quartile: ",quartile)
    # obtain the defined quartile of frequency per posID & grouping #
    sites <- arrange(sites, posID2, value, desc(freq))
    quartiles <- with(sites,
                      tapply(freq, posID2, quantile, probs=quartile, names=FALSE))
    sites$belowQuartile <- with(sites,freq < quartiles[posID2])
    rm(quartiles)
    
    if(any(sites$belowQuartile)) {
      # for values belowQuartile, see if any within defined windowSize of aboveQuartile #
      pos.be <- with(subset(sites,belowQuartile,drop=TRUE),
                     GRanges(IRanges(start=value,width=1), seqnames=posID2, freq=freq))
      pos.ab <- with(subset(sites,!belowQuartile,drop=TRUE),
                     GRanges(IRanges(start=value,width=1), seqnames=posID2, freq=freq))
      pos.overlap <- as.data.frame(as.matrix(findOverlaps(pos.be,pos.ab,maxgap=windowSize,
                                                          ignore.strand=TRUE)))
      
      # for overlapping values, merge them with the biggest respective aboveQuartile site #
      pos.overlap$freq <- values(pos.ab[pos.overlap[,"subjectHits"]])$freq
      maxs <- with(pos.overlap, tapply(freq, as.character(queryHits), max))
      pos.overlap$isMax <- with(pos.overlap, freq == maxs[as.character(queryHits)])
      rm(maxs)
      
      pos.overlap <- subset(pos.overlap,isMax, drop=TRUE)
      
      # if there are >1 biggest respective aboveQuartile site, then choose the closest one
      # if tied, then use the latter to represent the site
      counts <- xtabs(isMax~queryHits,pos.overlap)
      if(length(table(counts))>1) {
        toFix <- as.numeric(names(which(counts>1)))
        rows <- pos.overlap$queryHits %in% toFix
        pos.overlap$aboveQuartileValue <- 
          pos.overlap$belowQuartileValue <- pos.overlap$valueDiff <- 0
        pos.overlap$aboveQuartileValue[rows] <- 
          start(pos.ab[pos.overlap[rows,"subjectHits"]])
        pos.overlap$belowQuartileValue[rows] <- 
          start(pos.be[pos.overlap[rows,"queryHits"]])
        pos.overlap$valueDiff[rows] <- with(pos.overlap[rows,],
                                            abs(aboveQuartileValue-belowQuartileValue))
        mins <- with(pos.overlap[rows,], tapply(valueDiff, 
                                                as.character(queryHits), min))
        pos.overlap$isClosest <- TRUE
        pos.overlap$isClosest[rows] <- with(pos.overlap[rows,], 
                                            valueDiff == mins[as.character(queryHits)])
        pos.overlap <- subset(pos.overlap, isMax & isClosest,drop=TRUE)
        rm("counts","mins")
      }
      
      # trickle the changes back to the original dataframe#     
      pos.overlap$clusteredValue <- start(pos.ab[pos.overlap[,"subjectHits"]])    
      pos.overlap$posID2 <- as.character(seqnames(pos.be[pos.overlap[,"queryHits"]]))
      
      # for cases where no overlap was found, try clustering to themselves #
      rows <- which(!1:length(pos.be) %in% pos.overlap$query)
      loners <- pos.be[rows]
      if(length(loners)>0) {
        times.rep <- values(loners)[["freq"]]
        res <- clusterSites(rep(as.character(seqnames(loners)), times=times.rep),
                            rep(start(loners),times=times.rep),
                            byQuartile=FALSE)
      }
      pos.overlap <- rbind(pos.overlap[,c("queryHits","clusteredValue")], 
                           data.frame(queryHits=rows, 
                                      clusteredValue=as.numeric(res$clusteredValue)))
      
      sites$clusteredValue <- sites$value
      sites$clusteredValue[sites$belowQuartile][pos.overlap[,"queryHits"]] <- 
        pos.overlap$clusteredValue
      stopifnot(any(!is.na(sites$clusteredValue)))            
    } else {
      message("No sites found below defined quartile. Try to increase ",
              "the quartile or use standard clustering, byQuartile=FALSE.")
    }
  } else {
    message("Clustering by minimum overlap.")
    
    sites <- split(sites, sites$grouping)
    
    sites <- foreach(x=iter(sites), .inorder=FALSE, 
                     .packages=c("GenomicRanges","plyr"), 
                     .combine=rbind) %dopar% {
      
      ## find overlapping positions using findOverlaps() using maxgap adjusted by windowSize! ##
      sites.gr <- with(x, GRanges(seqnames=posID2, IRanges(start=value, width=1), 
                                  strand="*", freq))
      
      # the key part is ignoreSelf=TRUE,ignoreRedundant=FALSE...helps overwrite values at later step
      res <- as.data.frame(as.matrix(findOverlaps(sites.gr, ignoreSelf=TRUE, 
                                                  ignoreRedundant=FALSE,
                                                  select="all", maxgap=windowSize))) 
      if(nrow(res)>0) {
        # add accessory columns to dictate decision making!
        # q = query, s = subject, val = value, freq = frequency of query/subject
        res$q.val <- start(sites.gr)[res$queryHits]
        res$s.val <- start(sites.gr)[res$subjectHits]
        res$q.freq <- sites.gr$freq[res$queryHits]
        res$s.freq <- sites.gr$freq[res$subjectHits]
        res$dist <- with(res,abs(q.val-s.val))
        stopifnot(!any(res$dist>windowSize)) ## do safety checking!
        
        # favor a lower value where frequence/cloneCount is tied, else use the value of the highest frequency!
        res$val <- with(res,ifelse(q.freq==s.freq, 
                                   ifelse(q.val < s.val,q.val,s.val), 
                                   ifelse(q.freq >= s.freq,q.val,s.val))) 
        
        # For cases where there are >1 matches between query & subject...
        # find the one with the highest frequency and merge with that.
        # If all frequencies are the same, then use the lowest value to represent the cluster!
        res$maxFreq <- with(res, pmax(q.freq, s.freq))    
        maxes <- with(res, tapply(maxFreq, queryHits, max))
        res$ismaxFreq <- with(res, maxFreq==maxes[as.character(queryHits)])        
        
        ## VIP step...this is what merges high value to low value for ties in the hash structure below!!!
        res <- arrange(res, desc(queryHits), val)
        clustered <- unique(subset(res,ismaxFreq)[,c("queryHits","val")])
        clustered <- with(clustered, split(val, queryHits))
        
        ## make sure there is only one entry per hit
        ## this is useful in situations when multiple query & subject are off by 1bp
        ## i.e. queryHits: 1,2,3,4; subjectHits: 1,2,3,4; vals: 31895692 31895693 31895694 31895695
        clustered <- unlist(sapply(clustered, "[[", 1))        
        
        # trickle results back to sites
        x$clusteredValue <- x$value
        x$clusteredValue[as.numeric(names(clustered))] <- as.numeric(clustered)
        rm("clustered","res")
        cleanit <- gc()
      } else {
        message("No locations found within ", windowSize, 
                "bps for ",x$grouping[1],"...no clustering performed!")
        x$clusteredValue <- x$value
      }    
      x
    }        
  }
  
  message("\t - Adding clustered value frequencies.")
  # get frequency of clusteredValue
  counts <- count(sites[,-grep("value",names(sites),fixed=TRUE)],
                  c("posID2","clusteredValue"), wt_var="freq")
  names(counts)[grep("freq",names(counts),fixed=TRUE)] <- "clusteredValue.freq"
  sites <- merge(sites,counts)
  
  if(byQuartile) {
    sites <- sites[,c("posID","value","freq","clusteredValue",
                      "clusteredValue.freq","grouping")]
  }
  
  sites$posID2<-NULL
  if(is.null(grouping)) { sites$grouping<-NULL }
  if(is.null(weight)) { sites$weight<-NULL }
  
  return(sites)
}

#' Make OTUs of discrete positions grouped by reads
#'
#' Given a group of discrete positions per read/clone, the function tries to yield a unique OTU (operation taxinomical unit) ID for the collection based on overlap of discrete positions to other reads/clones by grouping. This is mainly useful when each readID has many posID which needs to be considered as one single group of sites.
#'
#' @param posID a vector of discrete positions, i.e. Chr,strand,Position.
#' @param readID a vector of read/clone names which is unique to each row, i.e. deflines.
#' @param grouping additional vector of grouping by which to pool the rows (i.e. samplenames). Default is NULL.
#' @param psl.rd a RangedData/GRanges object returned from \code{\link{clusterSites}}. Default is NULL. 
#'
#' @note The algorithm for making OTUs of sites is as follows: for each readID check how many positions are there. Separate readIDs with only position from the rest. Check if any readIDs with >1 position match to any readIDs with only one position. If there is a match, then assign both readIDs with the same OTU ID. Check if any positions from readIDs with >1 position match any other readIDs with >1 position. If yes, then assign same OTU ID to all readIDs sharing 1 or more positions. 
#'
#' @return a data frame with posID, readID, grouping, and otuID. If psl.rd parameter is defined, then a RangedData/GRanges object where object is first filtered by clusterTopHit column and the otuID column appended at the end.
#'
#' @seealso \code{\link{clusterSites}}, \code{\link{otuSites2}}, \code{\link{crossOverCheck}}, \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #otuSites(posID=c('chr1-1000','chr1-1000','chr2-1000','chr2+1000','chr15-1000','chr16-1000','chr11-1000'), readID=paste('read',sample(letters,7),sep='-'), grouping=c('a','a','a','b','b','b','c'))
#' #otuSites(psl.rd=test.psl.rd)
#'
otuSites <- function(posID=NULL, readID=NULL, grouping=NULL, 
                     psl.rd=NULL, parallel=TRUE) {
  if(!parallel) { registerDoSEQ() }
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(readID))
  } else {
    
    isRangedData <- FALSE
    if(is(psl.rd, "RangedData")) {
      isRangedData <- TRUE
      psl.rd <- as(psl.rd, "GRanges")
    }
    
    ## find the otuID by clusters ##
    if(!"clusterTopHit" %in% colnames(mcols(psl.rd)) | 
         !"clusteredPosition" %in% colnames(mcols(psl.rd))) {
      stop("The object supplied in psl.rd parameter does not have 'clusterTopHit' ",
           "or 'clusteredPosition' column in it. Did you run clusterSites() on it?")
    }
    
    good.rows <- mcols(psl.rd)$clusterTopHit
    posIDs <- paste0(as.character(seqnames(psl.rd)), as.character(strand(psl.rd)), 
                     mcols(psl.rd)$clusteredPosition)
    if("qName" %in% colnames(mcols(psl.rd))) {
      readID <- mcols(psl.rd)$qName
    } else if ("Sequence" %in% colnames(mcols(psl.rd))) {
      readID <- mcols(psl.rd)$Sequence
    } else {
      stop("No readID type column found in psl.rd object.")
    }
    grouping <- if(is.null(grouping)) { rep("A",length(psl.rd)) } else { grouping }
    
    otus <- otuSites(posIDs[good.rows], readIDs[good.rows], grouping[good.rows])
    
    message("Adding otuIDs back to psl.rd.")        
    otuIDs <- with(otus, split(otuID, paste(posID, readID, grouping)))
    mcols(psl.rd)$otuIDs <- NA
    mcols(psl.rd)$otuIDs[good.rows] <- as.numeric(otuIDs[paste(posIDs, readIDs, 
                                                               grouping)[good.rows]])
    
    rm("otus","otuIDs","posIDs","readIDs","grouping")
    cleanit <- gc()
    
    if(isRangedData) {
      psl.rd <- as(psl.rd, "RangedData")
    }
    
    return(psl.rd)
  }
  
  groups <- if(is.null(grouping)) { "A" } else { grouping }
  sites <- data.frame(posID, readID, grouping=groups, stringsAsFactors=FALSE)
  rm(groups)
  
  ## get unique posIDs per readID by grouping 
  # use tapply instead of ddply() or by() because it's a lot faster on larger datasets #
  counts <- with(sites, tapply(posID, paste0(grouping,readID), 
                               function(x) {
                                 uniques <- sort(unique(x))
                                 list(paste(uniques,collapse=","), length(uniques))
                               }))
  reads <- unique(sites[,c("grouping","readID")])
  reads$posIDs <- sapply(counts[with(reads,paste0(grouping,readID))],"[[", 1)
  reads$counts <- sapply(counts[with(reads,paste0(grouping,readID))],"[[", 2)
  reads <- arrange(reads, grouping, posIDs)
  
  # create initial otuID by assigning a numeric ID to each collection of posIDs per grouping
  reads$otuID <- unlist(lapply(lapply(with(reads, split(posIDs,grouping)),
                                      as.factor),
                               as.numeric)) 
  reads$newotuID <- reads$otuID 
  reads$check <- TRUE
  
  ## see if readID with a unique or single posID matches up to a readID with >1 posIDs, if yes then merge
  singles <- reads$counts==1
  if(any(singles)) {
    message('Merging non-singletons with singletons.')
    toCheck <- with(reads[singles,], split(posIDs,grouping))
    toCheck.ids <- with(reads[singles,], split(otuID,grouping))
    toCheck <- sapply(names(toCheck), 
                      function(x) { 
                        names(toCheck[[x]]) <- toCheck.ids[[x]]; toCheck[[x]] 
                        })
    rm(toCheck.ids)
    
    allposIDs <- with(reads[!singles,],split(posIDs,grouping))
    
    for(f in intersect(names(toCheck), names(allposIDs))) {
      # this is crucial to avoid matching things like xyzABC to xyz
      query <- structure(paste0(toCheck[[f]],","), 
                         names=names(toCheck[[f]])) 
      subject <- paste0(allposIDs[[f]],",") 
      
      res <- sapply(query, grep, x=subject, fixed=TRUE)
      res <- res[sapply(res,length)>0]                    
      res <- structure(unlist(res,use.names=F),
                       names=rep(names(res), sapply(res,length)))
      reads[!singles & reads$grouping==f,"newotuID"][as.numeric(res)] <- 
        as.numeric(names(res))
      reads[reads$grouping==f,"check"][as.numeric(res)] <- FALSE
    }    
  }
  
  ## see if readIDs with >1 posID overlap with other readIDs of the same type ##
  ## this is useful when no readIDs were found with a unique or single posID ##
  ## merge OTUs with overlapping positions within same grouping ##
  message('Merging non-singletons.')
  reads$grouping <- as.character(reads$grouping)
  reads$posIDs <- as.character(reads$posIDs)
  reads <- split(reads, reads$grouping)
  reads <- foreach(x=iter(reads), .inorder=FALSE, .combine=rbind) %dopar% {		
    for(f in 1:nrow(x)) {
      if (x$check[f]) {
        posId <- x$posIDs[f]
        oldid <- x$otuID[f]
        tocheck <- x[-f,][x$check,]
        res <- which(unlist(lapply(lapply(strsplit(tocheck$posIDs,","), "%in%", 
                                          unlist(strsplit(posId,","))), any)))
        if(length(res)>0) {
          id <- tocheck$otuID[res]
          x[x$otuID %in% id,"check"] <- FALSE
          x[x$otuID %in% id,"newotuID"] <- oldid
        }
      }
    }
    x$check <- NULL
    x
  }
  
  ## trickle the OTU ids back to sites frame ##    
  ots.ids <- with(reads, split(newotuID, paste0(readID,grouping)))
  if(!is.numeric(ots.ids)) {
    stop("Something went wrong merging non-singletons. ",
         "Multiple OTUs assigned to one readID most likely!")
  }
  sites$otuID <- as.numeric(unlist(ots.ids[with(sites, paste0(readID,grouping))]))
  
  stopifnot(any(!is.na(sites$otuID)))
  rm(reads)
  cleanit <- gc()
  
  if(is.null(grouping)) { sites$grouping<-NULL }
  return(sites)
}

#' Bin values or make OTUs by assigning a unique ID to them within discrete factors.
#'
#' Given a group of values or genomic positions per read/clone, the function tries to yield a unique OTU (operation taxinomical unit) ID for the collection based on overlap of locations to other reads/clones by grouping. This is mainly useful when each read has many locations which needs to be considered as one single group of sites.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand). Required if psl.rd parameter is not defined.
#' @param value a vector of integer locations/positions that needs to be binned, i.e. genomic location. Required if psl.rd parameter is not defined. 
#' @param readID a vector of read/clone names which is unique to each row, i.e. deflines.
#' @param grouping additional vector of grouping by which to pool the rows (i.e. samplenames). Default is NULL.
#' @param psl.rd a RangedData/GRanges object returned from \code{\link{clusterSites}}. Default is NULL. 
#' @param maxgap max distance allowed between two non-overlapping position to trigger the merging. Default is 5.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}. Process is split by the grouping the column.
#'
#' @note The algorithm for making OTUs of sites is as follows: 
#' \itemize{
#'  \item for each grouping & posID, fix values off by maxgap parameter
#'  \item create bins of fixed values per readID
#'  \item assign arbitrary numeric ID to each distinct bins above & obtain its frequency
#'  \item perform overlap b/w readIDs with only one value (singletons) to readIDs with >1 value (non-singletons)
#'  \item   - for any overlapping values, tag non-singleton readID with the ID of singleton readID
#'  \item   - if non-singleton readID matched with more than one singleton readID, then pick on at random
#'  \item for any non-tagged & non-singleton readIDs, perform an overlap of values within themselves using the maxgap parameter
#'  \item   - tag any overlapping positions across any readID with the ID of most frequently occuring bin
#'  \item positions with no overlap are left as is with the original arbitrary ID
#' }
#' 
#' @return a data frame with binned values and otuID shown alongside the original input. If psl.rd parameter is defined, then a RangedData/GRanges object.
#'
#' @seealso \code{\link{clusterSites}}, \code{\link{otuSites}}, \code{\link{crossOverCheck}}, \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #otuSites2(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), value=c(1000,1003,5832,1000,12324,65738,928042), readID=paste('read',sample(letters,7),sep='-'), grouping=c('a','a','a','b','b','b','c'))
#' #otuSites2(psl.rd=test.psl.rd)
#'
otuSites2 <- function(posID=NULL, value=NULL, readID=NULL, grouping=NULL, 
                      psl.rd=NULL, maxgap=5, parallel=TRUE) {
  if(!parallel) { registerDoSEQ() }
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(value))
    stopifnot(!is.null(readID))
  } else {
    
    isRangedData <- FALSE
    if(is(psl.rd, "RangedData")) {
      isRangedData <- TRUE
      psl.rd <- as(psl.rd, "GRanges")
    }
    
    ## find the otuID by Positions ##
    isthere <- grepl("clusteredPosition", colnames(mcols(psl.rd)), ignore.case=TRUE)
    if(!any(isthere)) {
      message("No 'clusteredPosition' column found in the data. Using 'Position' as an alternative to use as values.")
      isthere <- grepl("Position", colnames(mcols(psl.rd)), ignore.case=TRUE)  
      if(!any(isthere)) {
        stop("The object supplied in psl.rd parameter does not have ",
             "'clusteredPosition' or 'Position' column in it. ",
             "Did you run getIntegrationSites() or clusterSites() on it?")
      }
    }

    good.rows <- TRUE
    value <- mcols(psl.rd)[[which(isthere)]]
    posID <- paste0(as.character(seqnames(psl.rd)), as.character(strand(psl.rd)))
    if("qName" %in% colnames(mcols(psl.rd))) {
      readID <- mcols(psl.rd)$qName
    } else if ("Sequence" %in% colnames(mcols(psl.rd))) {
      readID <- mcols(psl.rd)$Sequence
    } else {
      stop("No readID type column found in psl.rd object.")
    }
    
    grouping <- if(is.null(grouping)) { rep("A",length(psl.rd)) } else { grouping }
    
    otus <- otuSites2(posID=posID[good.rows], 
                      value=value[good.rows], 
                      readID=readID[good.rows], 
                      grouping=grouping[good.rows], 
                      parallel=parallel)
    
    message("Adding otuIDs back to psl.rd.")        
    otuIDs <- with(otus,
                   split(otuID,
                         paste0(posID,value,readID,grouping)))
    mcols(psl.rd)$otuID <- NA
    mcols(psl.rd)$otuID[good.rows] <- as.numeric(otuIDs[paste0(posID, value, readID, 
                                                               grouping)[good.rows]])
    
    message("Cleaning up!")
    rm("otus","otuIDs","value","posID","readID","grouping","good.rows")
    cleanit <- gc()
    
    if(isRangedData) {
      psl.rd <- as(psl.rd, "RangedData")
    }
    
    return(psl.rd)
  }
  
  groups <- if(is.null(grouping)) { "A" } else { grouping }
  
  ## fix values off by maxgap parameter for sanity! ##
  sites.clustered <- clusterSites(posID, value, groups, windowSize=maxgap, 
                                  parallel=parallel)
  
  sites <- data.frame(posID, value, readID,                      
                      grouping=groups, stringsAsFactors=FALSE)
  sites <- merge(sites, sites.clustered)
  sites$posID2 <- with(sites, paste0(posID,clusteredValue))
  rm("groups","sites.clustered")
  
  ## get unique positions per readID by grouping 
  # use tapply instead of ddply() or by() because it's a lot faster on larger datasets #
  counts <- with(sites, tapply(posID2, paste0(grouping,readID), 
                               function(x) {
                                 uniques <- sort(unique(x))
                                 list(paste(uniques,collapse=","), length(uniques))
                               }))
  reads <- unique(sites[,c("grouping","readID")])
  reads$posIDs <- sapply(counts[with(reads,paste0(grouping,readID))],"[[", 1)
  reads$counts <- sapply(counts[with(reads,paste0(grouping,readID))],"[[", 2)
  
  # create initial otuID by assigning a numeric ID to each collection of posIDs per grouping
  reads$otuID <- unlist(
    lapply(lapply(with(reads,split(posIDs,grouping)), as.factor), as.numeric)
  ) 
  sites <- merge(arrange(sites,grouping,readID), 
                 arrange(reads[,c("grouping","readID","counts","otuID")],
                         grouping,readID), 
                 by=c("grouping","readID"), all.x=TRUE)
  sites$posID2 <- NULL
  rm(reads)
  sites <- arrange(sites, grouping, posID, clusteredValue)
  sites$readID <- as.character(sites$readID)
  sites$grouping <- as.character(sites$grouping)
  
  sites.gr <- with(sites, GRanges(seqnames=posID, IRanges(start=clusteredValue,
                                                          width=1), 
                                  strand="*", readID, grouping, counts, 
                                  otuID, newotuID=otuID))
  mcols(sites.gr)$grouping <- as.character(mcols(sites.gr)$grouping)
  mcols(sites.gr)$readID <- as.character(mcols(sites.gr)$readID)
  mcols(sites.gr)$check <- TRUE
  sites.gr <- sort(sites.gr)
  cleanit <- gc()
  
  ## see if readID with a unique/single location matches up to a readID with >1 location, if yes then merge
  mcols(sites.gr)$singles <- mcols(sites.gr)$counts==1
  if(any(mcols(sites.gr)$singles)) {
    message('Merging non-singletons with singletons if any...')    
    sites.gr.list <- split(sites.gr, mcols(sites.gr)$grouping)
    sites.gr <- foreach(x=iter(sites.gr.list), .inorder=FALSE, .export="maxgap",
                        .packages="GenomicRanges", .combine=c) %dopar% {
                          sigs <- subset(x, mcols(x)$singles)
                          nonsigs <- subset(x, !mcols(x)$singles)
                          res <- findOverlaps(nonsigs, sigs, maxgap=maxgap, 
                                              select="first")
                          rows <- !is.na(res)
                          if(any(rows)) {
                            res <- data.frame(queryHits=which(rows), 
                                              subjectHits=res[rows])
                            res$sigsOTU <- mcols(sigs)$otuID[res$subjectHits]
                            
                            res$sigsReadID <- mcols(sigs)$readID[res$subjectHits]
                            res$nonsigsReadID <- mcols(nonsigs)$readID[res$queryHits]

                            res$sigPosID <- paste0(as.character(seqnames(sigs)),
                                                       start(sigs))[res$subjectHits]
                            res$nonsigPosID <- paste0(as.character(seqnames(nonsigs)),
                                                     start(nonsigs))[res$queryHits]
                            
                            ## if >1 OTU found per nonsigsReadID...choose lowest ID ## 
                            bore <- with(res, split(sigsOTU, nonsigsReadID))
                            bore <- sapply(sapply(sapply(bore, unique, simplify=FALSE), 
                                                  sort, simplify=FALSE), 
                                           "[[", 1)
                            res$OTU <- bore[res$nonsigsReadID]                                                        
                            
                            ## if >1 OTU found per nonsigPosID...choose lowest ID ##                          
                            res$OTU2 <- res$OTU
                            counts <- with(res, tapply(OTU2, nonsigPosID, 
                                                       function(x) length(unique(x))))                              
                            totest <- names(which(counts>1))
                            while(length(totest)>0) { 
                              #print(length(totest))
                              for(i in totest) {
                                rows <- res$nonsigPosID == i
                                rows <- res$nonsigsReadID %in% res$nonsigsReadID[rows]
                                rows <- rows | res$nonsigPosID %in% res$nonsigPosID[rows]
                                res$OTU2[rows] <- min(res[rows,"OTU"])
                              }
                              counts <- with(res, tapply(OTU2, nonsigPosID, 
                                                         function(x) length(unique(x))))                              
                              totest <- names(which(counts>1))
                            }

                            bore <- sapply(with(res, split(OTU2, nonsigsReadID)),
                                           unique)
                            if(!is.numeric(bore)) { 
                              ## safety check incase >1 OTU found per readID
                              bore <- sapply(bore, min)
                            }
                            rows <- mcols(nonsigs)$readID %in% names(bore)
                            mcols(nonsigs)[rows,"newotuID"] <- 
                              bore[mcols(nonsigs)[rows,"readID"]]
                            mcols(nonsigs)[rows,"check"] <- FALSE
                            mcols(sigs)[mcols(sigs)$readID %in% res$sigsReadID, "check"] <- FALSE
                          }
                          c(sigs,nonsigs)
                        }
    rm(sites.gr.list)
  }
  mcols(sites.gr)$singles <- NULL
  
  ## perform non-singletons overlap of values within maxgap ##
  ## merge OTUs with overlapping positions within same grouping ##
  message('Performing non-singletons overlap...')
  goods <- subset(sites.gr, !mcols(sites.gr)$check)
  sites.gr <- subset(sites.gr, mcols(sites.gr)$check)
  sites.gr.list <- split(sites.gr, mcols(sites.gr)$grouping)
  sites.gr <- foreach(x=iter(sites.gr.list), .inorder=FALSE, .export="maxgap",
                      .packages="GenomicRanges", .combine=c) %dopar% {		    
                        res <- findOverlaps(x, maxgap=maxgap, ignoreSelf=TRUE,
                                            ignoreRedundant=TRUE, select="all")
                        if(length(res)>0) {
                          res <- as.data.frame(res)
                          res$queryOTU <- mcols(x)$otuID[res$queryHits]
                          res$subjectOTU <- mcols(x)$otuID[res$subjectHits]
                          res$subjectReadID <- mcols(x)$readID[res$subjectHits]
                          res$subjectPosID <- paste0(as.character(seqnames(x)),
                                                     start(x))[res$subjectHits]
                          res$queryPosID <- paste0(as.character(seqnames(x)),
                                                   start(x))[res$queryHits]
                          
                          ## if >1 OTU found per subjectReadID...choose lowest ID ## 
                          bore <- with(res, split(queryOTU, subjectReadID))
                          bore <- sapply(sapply(sapply(bore, unique, simplify=FALSE), 
                                                sort, simplify=FALSE), 
                                         "[[", 1)
                          res$OTU <- bore[as.character(res$subjectReadID)]                                                    
                          
                          ## if >1 OTU found per subjectPosID...choose lowest ID ##                          
                          res$OTU2 <- res$OTU                        
                          counts <- with(res, tapply(OTU2, subjectPosID, 
                                                     function(x) length(unique(x))))                          
                          totest <- names(which(counts>1))
                          
                          while(length(totest)>0) {
                            #print(length(totest))
                            for(i in totest) {
                              rows <- res$subjectPosID == i 
                              rows <- rows | res$subjectReadID %in% res$subjectReadID[rows]                               
                              rows <- rows | res$subjectOTU %in% res$subjectOTU[rows]
                              rows <- rows | res$subjectPosID %in% res$subjectPosID[rows]
                              rows <- rows | res$queryOTU %in% res$subjectOTU[rows]
                              
                              res$OTU2[rows] <- min(res[rows,"OTU2"])
                            }                            
                            counts <- with(res, tapply(OTU2, subjectPosID, 
                                                       function(x) length(unique(x))))                          
                            totest <- names(which(counts>1))
                          }
                          
                          bore <- sapply(with(res, split(OTU2, subjectHits)),
                                         unique)
                          ## safety check incase >1 OTU found per subjectHits
                          stopifnot(is.numeric(bore)) 
                          
                          rows <- as.numeric(names(bore))
                          mcols(x)[rows,"newotuID"] <- as.numeric(bore)
                          mcols(x)[rows, "check"] <- FALSE
                          
                          bore <- sapply(with(res, split(OTU2, subjectReadID)),
                                         unique)
                          if(!is.numeric(bore)) { 
                            ## safety check incase >1 OTU found per readID
                            bore <- sapply(bore, min)
                          }
                          rows <- mcols(x)$readID %in% names(bore)
                          mcols(x)[rows,"newotuID"] <- bore[mcols(x)[rows,"readID"]]
                          mcols(x)[rows,"check"] <- FALSE
                        }
                        x
                      }
  sites.gr <- c(sites.gr, goods)
  rm("sites.gr.list","goods")
  cleanit <- gc()
  
  ## trickle the OTU ids back to sites frame ##    
  ots.ids <- sapply(split(mcols(sites.gr)$newotuID,
                          paste0(mcols(sites.gr)$readID, 
                                 mcols(sites.gr)$grouping)),
                    unique)
  
  if(!is.numeric(ots.ids) | any(sapply(ots.ids, length)>1)) {
    stop("Something went wrong merging non-singletons. ",
         "Multiple OTUs assigned to one readID most likely!")
  }
  sites$otuID <- as.numeric(unlist(ots.ids[with(sites,
                                                paste0(readID,grouping))],
                                   use.names=F))
  
  stopifnot(any(!is.na(sites$otuID)))
  cleanit <- gc()
  
  if(is.null(grouping)) { sites$grouping <- NULL }
  
  sites$clusteredValue <- NULL; sites$clusteredValue.freq <- NULL
  sites$freq <- NULL; sites$counts <- NULL
  
  return(sites)
}

#' Bin values or make ISUs by assigning a unique ID to them within discrete factors.
#'
#' Given a group of values or genomic positions per read/clone, the function tries to yield a unique ISU (Integration Site Unit) ID for the collection based on overlap of locations to other reads/clones by grouping. This is mainly useful when each read has many locations which needs to be considered as one single group of sites.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand). Required if psl.rd parameter is not defined.
#' @param value a vector of integer locations/positions that needs to be binned, i.e. genomic location. Required if psl.rd parameter is not defined. 
#' @param readID a vector of read/clone names which is unique to each row, i.e. deflines.
#' @param grouping additional vector of grouping by which to pool the rows (i.e. samplenames). Default is NULL.
#' @param psl.rd a RangedData/GRanges object returned from \code{\link{clusterSites}}. Default is NULL. 
#' @param maxgap max distance allowed between two non-overlapping position to trigger the merging. Default is 5.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}. Process is split by the grouping the column.
#'
#' @note The algorithm for making isus of sites is as follows: for each readID check how many positions are there. Separate readIDs with only position from the rest. Check if any readIDs with >1 position match to any readIDs with only one position. If there is a match, then assign both readIDs with the same ISU ID. Check if any positions from readIDs with >1 position match any other readIDs with >1 position. If yes, then assign same ISU ID to all readIDs sharing 1 or more positions.
#'
#' @return a data frame with binned values and isuID shown alongside the original input. If psl.rd parameter is defined, then a RangedData/GRanges object where object is first filtered by clusterTopHit column and the isuID column appended at the end.
#'
#' @seealso \code{\link{clusterSites}}, \code{\link{isuSites}}, \code{\link{crossOverCheck}}, \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #isuSites2(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), value=c(rep(1000,2),5832,1000,12324,65738,928042), readID=paste('read',sample(letters,7),sep='-'), grouping=c('a','a','a','b','b','b','c'))
#' #isuSites2(psl.rd=test.psl.rd)
#'
isuSites <- function(posID=NULL, value=NULL, readID=NULL, grouping=NULL, 
                      psl.rd=NULL, maxgap=5, parallel=TRUE) {
  
  res <- otuSites2(posID=posID, value=value, readID=readID, grouping=grouping, 
                   psl.rd=psl.rd, maxgap=maxgap, parallel=parallel)
  
  if(is(res,"GRanges")) {
    cols <- grep('otu',colnames(mcols(res)))
    colnames(mcols(res))[cols] <- gsub("otu","isu",colnames(mcols(res))[cols])
  } else {
    cols <- grep('otu',colnames(res))
    colnames(res)[cols] <- gsub("otu","isu",colnames(res)[cols])
  }
  
  res
}

#' Remove values/positions which are overlapping between discrete groups based on their frequency.
#'
#' Given a group of discrete factors (i.e. position ids) and integer values, the function tests if they overlap between groups. If overlap is found, then the group having highest frequency of a given position wins, else the position is filtered out from all the groups. The main use of this function is to remove crossover sites from different samples in the data.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand). Required if psl.rd parameter is not defined.
#' @param value a vector of integer locations/positions that needs to be binned, i.e. genomic location. Required if psl.rd parameter is not defined. 
#' @param grouping additional vector of grouping by which to pool the rows (i.e. samplenames). Default is NULL.
#' @param weight a numeric vector of weights to use when calculating frequency of value by posID and grouping if specified. Default is NULL.
#' @param psl.rd a RangedData/GRanges object. Default is NULL. 
#'
#' @return a data frame of the original input with columns denoting whether a given row is crossover or not. If psl.rd parameter is defined, then a RangedData/GRanges object with 'isCrossover' column appended at the end.
#'
#' @seealso  \code{\link{clusterSites}}, \code{\link{otuSites}}, \code{\link{otuSites2}}, \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #crossOverCheck(posID=c('chr1-','chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), value=c(rep(1000,3),5832,1000,12324,65738,928042), grouping=c('a','a','b','b','b','b','c','c'))
#' #crossOverCheck(psl.rd=test.psl.rd)
#'
crossOverCheck <- function(posID=NULL, value=NULL, grouping=NULL, 
                           weight=NULL, psl.rd=NULL) {
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(value))
  } else {     
    
    isRangedData <- FALSE
    if(is(psl.rd, "RangedData")) {
      isRangedData <- TRUE
      psl.rd <- as(psl.rd, "GRanges")
    }
    
    if(!"clusterTopHit" %in% colnames(mcols(psl.rd)) | 
         !"clusteredPosition" %in% colnames(mcols(psl.rd))) {
      stop("The object supplied in psl.rd parameter does not have 'clusterTopHit' ",
           "or 'clusteredPosition' column in it. Did you run clusterSites() on it?")
    }
    
    posID <- paste0(as.character(seqnames(psl.rd)), as.character(strand(psl.rd)))
    if("clusteredPosition" %in% colnames(mcols(psl.rd))) {
      message("Using clusteredPosition column from psl.rd as the value parameter.")
      value <- mcols(psl.rd)$clusteredPosition
      good.row <- mcols(psl.rd)$clusterTopHit
    } else if("Position" %in% colnames(mcols(psl.rd))) {
      message("Using Position column from psl.rd as the value parameter.")
      value <- mcols(psl.rd)$Position
      good.row <- rep(TRUE, length(value))
    } else {
      message("Using start(psl.rd) as the value parameter.")
      value <- start(psl.rd)
      good.row <- rep(TRUE, length(value))
    }
    
    isthere <- grepl("isMultiHit", colnames(mcols(psl.rd)), ignore.case=TRUE)
    if(any(isthere)) { ## see if multihit column exists
      message("Found 'isMultiHit' column in the data. ",
              "These rows will be ignored for the calculation.")
      good.row <- good.row & !mcols(psl.rd)[[which(isthere)]]
    }
    
    if(is.null(weight)) { ## see if clonecount column exists
      isthere <- grepl("clonecount", colnames(mcols(psl.rd)))
      if(any(isthere)) { weight <- mcols(psl.rd)[[which(isthere)]] } else { weight <- weight }
    }
    grouping <- if(is.null(grouping)) { "" } else { grouping }
    crossed <- crossOverCheck(posID[good.row], value[good.row], 
                              grouping=grouping[good.row], weight=weight[good.row])
    
    message("Adding isCrossover data back to psl.rd.")        
    crossedValues <- with(crossed, split(isCrossover, paste0(posID,value,grouping)))
    if(any(sapply(crossedValues, length)>1)) {
      stop("Error in crossOverCheck: sampling culprits... ", 
           paste(names(crossedValues[which(sapply(crossedValues, length)>1)]), 
                 collapse=", "))
    }
    mcols(psl.rd)$isCrossover <- FALSE
    mcols(psl.rd)$isCrossover[good.row] <- 
      as.logical(crossedValues[paste0(posID,value,grouping)[good.row]])
    
    message("Cleaning up!")
    rm("posID","value","grouping","crossed","crossedValues")
    cleanit <- gc()
    
    if(isRangedData) {
      psl.rd <- as(psl.rd, "RangedData")
    }
    
    return(psl.rd)
  }
  
  # get frequencies of each posID & value combination by grouping #
  groups <- if(is.null(grouping)) { "" } else { grouping }
  weight2 <- if(is.null(weight)) { 1 } else { weight }
  sites <- data.frame(posID, value, grouping=groups, weight=weight2, 
                      stringsAsFactors=FALSE)
  sites <- arrange(sites, grouping, posID, value)
  sites <- count(sites, c("posID","value","grouping"), wt_var="weight")
  rm("groups","weight2")
  
  sites$isCrossover <- FALSE
  sites$FoundIn <- sites$grouping
  
  # find overlapping positions & pick the winner based on frequencies #
  sites.gr <- with(sites, GRanges(seqnames=posID, IRanges(start=value,width=1), 
                                  strand="*", grouping, freq))
  res <- findOverlaps(sites.gr, maxgap=1, ignoreSelf=TRUE, 
                      ignoreRedundant=FALSE, select="all")
  if(length(res)>0) {
    res <- as.data.frame(res)
    res$qgroup <- mcols(sites.gr)$grouping[res$queryHits]
    res$sgroup <- mcols(sites.gr)$grouping[res$subjectHits]
    res$qfreq <- mcols(sites.gr)$freq[res$queryHits]
    res$sfreq <- mcols(sites.gr)$freq[res$subjectHits]
    res <- ddply(res, .(queryHits), summarize, 
                 FoundIn=paste(sort(unique(c(qgroup,sgroup))),collapse=","),
                 isBest=all(qfreq>sfreq))
    sites$isCrossover[subset(res,!isBest)$queryHits] <- TRUE  
    sites$FoundIn[res$queryHits] <- res$FoundIn
  }
  sites
}

#' Find the integration sites and add results to SampleInfo object. 
#'
#' Given a SampleInfo object, the function finds integration sites for each sample using their respective settings and adds the results back to the object. This is an all-in-one function which aligns, finds best hit per read per sample, cluster sites, and assign ISU IDs. Depending on the aligner chosen, it calls \code{\link{blatSeqs}}, \code{\link{read.psl}}, \code{\link{getIntegrationSites}}, \code{\link{clusterSites}}, \code{\link{otuSites}}. There must be linkered reads within the sampleInfo object in order to use this function using the default parameters. If you are planning on BLATing non-linkered reads, then change the seqType to one of accepted options for the 'feature' parameter of \code{\link{extractSeqs}}, except for '!' based features.
#'
#' @param sampleInfo sample information SimpleList object outputted from \code{\link{findLinkers}}, which holds decoded, primed, LTRed, and Linkered sequences for samples per sector/quadrant along with metadata.
#' @param seqType which type of sequence to align and find integration sites for. Default is NULL and determined automatically based on type of restriction enzyme or isolation method used. If restriction enzyme is Fragmentase or Mu, then this parameter is set to genomicLinkered, else it is genomic. Any one of following options are valid: genomic, genomicLinkered, decoded, primed, LTRed, linkered.
#' @param aligner which aligner to use: BLAT(default) or Subread. Use 'Subread' for illumina data.
#' @param port a port number to host the gfServer with. Only used if aligner='BLAT'. Default is 5560.
#' @param host name of the machine running gfServer. Only used if aligner='BLAT'. Default is 'localhost'.
#' @param genomeIndices an associative character vector of freeze to full or relative path of respective of indexed genomes from BLAT(.nib or .2bit files) or Subreads(.subread.index). For example: c("hg18"="/usr/local/blatSuite34/hg18.2bit", "mm8"="/usr/local/blatSuite34/mm8.2bit"). Be sure to supply an index per freeze supplied in the sampleInfo object. Default is NULL.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param samplenames a vector of samplenames to process. Default is NULL, which processes all samples from sampleInfo object.
#'
#' @return a SimpleList object similar to sampleInfo parameter supplied with new data added under each sector and sample. New data attributes include: psl, and sites. The psl attributes holds the genomic hits per read along with QC information. The sites attribute holds the condensed integration sites where genomic hits have been clustered by the Position column and cherry picked to have each site pass all the QC steps. 
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#'
#' @seealso \code{\link{findPrimers}}, \code{\link{findLTRs}}, \code{\link{findLinkers}}, \code{\link{startgfServer}}, \code{\link{read.psl}}, \code{\link{blatSeqs}}, \code{\link{blatListedSet}}, \code{\link{pslToRangedObject}}, \code{\link{clusterSites}}, \code{\link{otuSites2}}, \code{\link{crossOverCheck}}, \code{\link{getIntegrationSites}}
#'
#' @export
#'
#' @examples 
#' #findIntegrations(sampleInfo)
#'
findIntegrations <- function(sampleInfo, seqType=NULL, aligner="BLAT",
                             port=5560, host="localhost", genomeIndices=NULL,
                             parallel=TRUE, samplenames=NULL) {    

  .checkArgs_SEQed()
  
  ## test if there are linkered sequences in the sampleinfo object if specific feature/seqType is not defined ##   
  feature <- ifelse(is.null(seqType), "linkered", seqType)
  message("Checking for ", feature, " reads.")		
  featured <- extractFeature(sampleInfo, feature=feature)
  samplesfeatured <- sapply(featured, names, simplify=FALSE)
  sectorsfeatured <- names(which(sapply(sapply(featured,length),">",0)))
  rm(featured)
  cleanit <- gc()
  
  if(length(sectorsfeatured)==0) {
    stop("No ", feature, " reads found in sampleInfo object provided.")
  }
  
  ## subset specific samples if defined ##
  samplesToProcess <- unlist(samplesfeatured,use.names=F)
  if(!is.null(samplenames)) {
    samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
  }
  
  message("Creating hashes of settings for aligning and processing.")
  
  ## setup settings hashes for aligning the genomic seqs by species & processing hits later ##
  for(setting in c("restrictionenzyme", "freeze", "startwithin", 
                   "alignratiothreshold", "genomicpercentidentity", 
                   "clustersiteswithin", "keepmultihits")) {
    setter <- extractFeature(sampleInfo, samplename=samplesToProcess, feature=setting)
    names(setter) <- NULL
    setter <- unlist(setter)
    assign(setting,setter)
  }
  
  ## Align by respective species ##
  pslFiles <- c()        
  for (f in unique(freeze)) {
    message("Aligning to: ", f)
    
    message("Getting sequences to align")        
    # get sequences to align #
    if (is.null(seqType)) {
      wanted <- names(restrictionenzyme[!grepl("FRAG|SONIC|MU",
                                               restrictionenzyme,ignore.case=TRUE)])
      wanted <- wanted[wanted %in% names(freeze[freeze==f])]
      seqs <- extractSeqs(sampleInfo, samplename=wanted, 
                          feature="genomic", minReadLength=5)
      if(any(as.numeric(sapply(seqs,length))>0)) {
        write.listedDNAStringSet(seqs, filePrefix=paste0("processed",f))
      }
      
      wanted <- names(restrictionenzyme[grepl("FRAG|SONIC|MU",
                                              restrictionenzyme,ignore.case=TRUE)])
      wanted <- wanted[wanted %in% names(freeze[freeze==f])]
      seqs <- extractSeqs(sampleInfo, samplename=wanted, 
                          feature="genomicLinkered", minReadLength=5)
      if(any(as.numeric(sapply(seqs,length))>0)) {
        write.listedDNAStringSet(seqs, filePrefix=paste0("processed",f))
      }
    } else {
      seqs <- extractSeqs(sampleInfo, samplename=names(freeze[freeze==f]), 
                          feature=seqType, minReadLength=5)
      if(any(as.numeric(sapply(seqs,length))>0)) {
        write.listedDNAStringSet(seqs, filePrefix=paste0("processed",f))
      }
    }
    
    # Align seqs #
    pslFile <- blatSeqs(query=paste0("processed",f,".*.fa$"), 
                        subject=genomeIndices[[f]], standaloneBlat=FALSE, 
                        host=host, port=port, parallel=parallel, gzipResults=TRUE)                
    
    message("Cleaning!")
    # add pslFiles for later use #
    pslFiles <- c(pslFiles,pslFile)
    cleanit <- gc()
    system(paste("rm",paste0("processed",f,".*.fa")))
  }
  
  message("Reading PSL files.")
  ## read all hits and split by samples ##
  psl <- read.psl(pslFiles, bestScoring=TRUE, asGRanges=TRUE, 
                  removeFile=TRUE, parallel=FALSE)
  cleanit <- gc()
  psl$setname <- sub("^(.+)-(.+)$","\\1", psl$qName)
  psl <- split(psl, psl$setname)
  
  ## begin processing hits ##
  psl.hits <- foreach(x=iter(names(psl)), .inorder=TRUE,
                      .export=c("psl", "sampleInfo", "startwithin", 
                                "alignratiothreshold", "genomicpercentidentity", 
                                "clustersiteswithin", "keepmultihits",
                                "getIntegrationSites", "clusterSites", 
                                "isuSites"),
                      .packages=c("GenomicRanges","plyr")) %dopar% {
    message("Processing ",x)
    
    # add qc info for bestscoring hits #
    psl.x <- getIntegrationSites(psl[[x]], startWithin=startwithin[[x]], 
                                 alignRatioThreshold=alignratiothreshold[[x]], 
                                 genomicPercentIdentity=genomicpercentidentity[[x]])
    
    # filter multihits if applicable #
    if(!as.logical(keepmultihits[[x]])) {
      psl.x <- psl.x[!psl.x$isMultiHit, ]
    }
    
    # cluster sites by positions #
    psl.x <- clusterSites(psl.rd=psl.x, windowSize=clustersiteswithin[[x]])
    
    # get sites ISU for tagging multihits #
    if(as.logical(keepmultihits[[x]])) {        
      psl.x <-isuSites(psl.rd=psl.x)
    }
    
    psl.x
  }
  names(psl.hits) <- names(psl)
  
  message("Adding PSL hits back to the object.")
  sampleInfo <- addFeature(sampleInfo, sector=NULL, samplename=names(psl.hits),
                           feature="psl", value=psl.hits)
  
  message("Adding sites back to the object.")
  psl.hits <- sapply(psl.hits, function(x) {
    x <- DataFrame(x)      
    x <- subset(x, clusterTopHit & pass.allQC, 
                select=setdiff(colnames(x), c('matches', 'misMatches', 'repMatches', 
                                              'nCount', 'qNumInsert', 'qBaseInsert', 
                                              'tNumInsert', 'tBaseInsert', 'tSize', 
                                              'blockCount', 'blockSizes', 'qStarts', 
                                              'score', 'tStarts', 'pass.startWithin',
                                              'alignRatio', 'pass.alignRatio', 
                                              'percIdentity', 'pass.percIdentity', 
                                              'pass.allQC', 'clusterTopHit', 'width',
                                              'Position')))
    x$start <- x$end <- x$clusteredPosition
    x$clusteredPosition <- NULL    
    GRanges(x)
  })
  
  sampleInfo <- addFeature(sampleInfo, sector=NULL, samplename=names(psl.hits),
                           feature="sites", value=psl.hits)
  
  cleanit <- gc()
  
  sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
  return(sampleInfo)
}

#' Simple summary of a sampleInfo object.
#'
#' Give a simple summary of major attributes in sampleInfo/SimpleList object.
#'
#' @param sampleInfo sample information SimpleList object, which samples per sector/quadrant information along with other metadata.
#'
#' @return a dataframe summarizing counts of major attributes per sample and sector. 
#'
#' @export
#'
#' @examples 
#' #summary.simple(sampleInfo)
#'
summary.simple <- function(sampleInfo) {
  stopifnot(is(sampleInfo,"SimpleList"))
  message("Total sectors:", paste(names(sampleInfo$sectors),collapse=","), "\n")
  do.call(rbind, lapply(names(sampleInfo$sectors), function(sector) {
    res.df <- data.frame(Sector=sector, 
                         SampleName=as.character(extractFeature(sampleInfo,
                                                                sector=sector,
                                                                feature="samplename")[[sector]]))
    res.df$SampleName <- as.character(res.df$SampleName)
    for (metaD in c("decoded","primed","LTRed","linkered","psl","sites")) {
      res <- extractFeature(sampleInfo, sector=sector, feature=metaD)[[sector]]
      
      if(is(res,"DataFrame")) {
        res <- t(sapply(as.list(res), sapply, length))
        if(length(res)>0) {
          res.df[,metaD] <- res[res.df$SampleName,]
        } else {
          res.df[,metaD] <- NA
        }
      } else {
        res <- sapply(res, length)
        if(length(res)>0) {
          res.df[,metaD] <- res[res.df$SampleName]
        } else {
          res.df[,metaD] <- NA
        }
      }
    }
    res.df        
  }))    
}

#' Elegant summary of a sampleInfo object.
#'
#' Give an elegant summary of all the attributes in sampleInfo/SimpleList object.
#'
#' @param sampleInfo sample information SimpleList object, which samples per sector/quadrant information along with other metadata.
#' @param samplenames a vector of samplenames to process. Default is NULL, which processes all samples from sampleInfo object.
#'
#' @export
#'
#' @examples 
#' #summary.elegant(sampleInfo)
#'
summary.elegant <- function(sampleInfo, samplenames=NULL) {
  stopifnot(is(sampleInfo,"SimpleList"))
  cat("Names of all sectors processed:",
      paste(names(sampleInfo$sectors),collapse=","), 
      "\n")
  for(sector in names(sampleInfo$sectors)) {
    cat(rep("*",50),"\n")
    cat("Sector:",sector,"\n")
    
    sampleList <- as.character(extractFeature(sampleInfo,sector=sector,
                                              feature="samplename")[[sector]])
    isPaired <- extractFeature(sampleInfo, sector=sector, 
                               feature="pairedend")[[sector]]
    if(!is.null(samplenames)) {
      sampleList <- sampleList[sampleList %in% samplenames]
    }
    
    for(sample.i in sampleList) {
      cat("\t",sample.i,":\n\t")
      bore <- strsplit(extractFeature(sampleInfo,
                                      sector=sector,
                                      samplename=sample.i,
                                      feature="metadata")[[sector]], ", ")
      metadatalist <- as.character(unlist(bore))
      totalDecoded <- extractFeature(sampleInfo, sector=sector, 
                                     samplename=sample.i, feature="decoded")[[sector]]
      
      decode.test <- if(isPaired[[sample.i]]) {
        all(sapply(totalDecoded, sapply, length)>0)
      } else {
        length(totalDecoded)>0
      }
      
      if(decode.test) {
        totalDecoded <- if(isPaired[[sample.i]]) {
          sapply(totalDecoded, sapply, length) 
        } else { 
          length(totalDecoded[[1]])
        }
        
        meta.processed <- grep("ed$", metadatalist, value=TRUE)
        meta.sites <- c("psl", "sites")
        for(metadata.i in metadatalist) {
          cat("| ", metadata.i, ": ")
          feature.i <- extractFeature(sampleInfo, sector=sector, samplename=sample.i,
                                      feature=metadata.i)[[sector]][[1]]
          
          counts <- if(is.list(feature.i)) {
            sapply(feature.i, length)
          } else {
            length(feature.i)
          }
          
          if(metadata.i %in% meta.processed) {
            cat(length(feature.i),"(",
                round(100*(counts/totalDecoded),2),"%) ",
                sep="","\n")
          } else if(metadata.i %in% meta.sites) {
            counts <- if(is.list(feature.i)) {
              sapply(feature.i, function(x) length(unique(x$qName)))
            } else {
              length(unique(feature.i$qName))
            }
            cat(counts,"(",
                round(100*(counts/totalDecoded),2),"%) ",
                sep="","\n")
          } else {
            cat(feature.i," ")
          }                
        }
      } else {
        cat("\t None decoded. \n")
      }
      cat("\n",rep("~",40),"\n")
    }
  }    
}
