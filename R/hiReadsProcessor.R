#' Functions for process LM-PCR reads with hiReadsProcessor.
#'
#' hiReadsProcessor contains set of functions which allow users to process single-end LM-PCR sequence data coming out of the 454/Solexa sequencer. Given an excel file containing parameters for demultiplexing and sample metadata, the functions automate trimming of adaptors and identification of the genomic product. In addition, if IntSites MySQL database is setup, the sequence attrition is loaded into respective tables for post processing setup and analysis.
#'
#' @import Biostrings GenomicRanges parallel doBy RMySQL xlsx plyr dataframe
#' @docType package
#' @name hiReadsProcessor
NULL

#' Read contents of a sequencing folder and make a SimpleList object
#'
#' Given a sequencing folder path, sample information file path, and sequence file extension pattern, the function returns a list of variables required to process the data. The function also calls \code{\link{read.sampleInfo}} which reads in sample processing metadata and formats it if needed.
#'
#' @param SequencingFolderPath full or relative path to the sequencing folder
#' @param sampleInfoFilePath full or relative path to the sample information file, which holds samples to quadrant/lane associations along with other metadata required to trim sequences or process it. Default to NULL, where the function tries to find xls or tab deliminated txt file in the sequencing folder which sounds similar to 'sampleinfo' and present you with choices of file to select from.
#' @param seqfilePattern regex to describe sequence file endings. See examples. Default is "\\.TCA.454Reads.fna$".
#' @param interactive whether to prompt each time the function encounters an issue, or use the defaults. Default is TRUE.
#'
#' @return a SimpleList list which is used by other functions to process and decode the data.
#'
#' @note one must make sure that each sequencing file has sector name/number prefixed at the beginning, else \code{\link{decodeByBarcode}} will fail trying to find the filename.
#'
#' @seealso \code{\link{read.sampleInfo}}, \code{\link{decodeByBarcode}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples  
#' #read.SeqFolder("~/Downloads/454Runs/2011_07_22",seqfilePattern=".+fna$")
#' #read.SeqFolder(".",seqfilePattern=".+fastq$")
#'
read.SeqFolder <- function(SequencingFolderPath=NULL, sampleInfoFilePath=NULL, seqfilePattern="\\.TCA.454Reads.fna$", interactive=TRUE) {
    if(is.null(SequencingFolderPath)) {
        stop("No Sequencing Folder Path provided.")
    }
    
    ## get the SequencingFolderPath right!
    SequencingFolderPath <-  normalizePath(SequencingFolderPath,mustWork=TRUE)
    
    seqFilePaths <- list.files(path=SequencingFolderPath, recursive=TRUE, full.names=TRUE, pattern=seqfilePattern)
    if(length(seqFilePaths)==0) {
        stop(paste("No files found in the folder",SequencingFolderPath,"matching following pattern:",seqfilePattern))
    }
    
    ## read the sample info file
    if(is.null(sampleInfoFilePath)) {
        possibleFiles <- list.files(path=SequencingFolderPath, recursive=TRUE, full.names=TRUE, pattern=".*sampleinfo.+", ignore.case=TRUE)
        if (length(possibleFiles)==0) {
            stop(paste("No sample information file found in folder :",SequencingFolderPath))
        } else {
            if(interactive & length(possibleFiles)>1) {
                message("Please choose a sample information file to read the meta data from:\n",paste(1:length(possibleFiles),possibleFiles,sep=": ",collapse="\n"))
                choice <- scan(what=integer(0),n=1,quiet=TRUE,multi.line=FALSE)
            } else {
                choice <- 1            
            }
            message("Choosing ",possibleFiles[choice]," as sample information file.")
        }
        sampleInfoFilePath <- possibleFiles[choice]
    }
    sampleInfo <- read.sampleInfo(sampleInfoFilePath,interactive=interactive)
    
    if(length(sampleInfo)!=length(seqFilePaths)) {
        warning("Number of sectors (",length(sampleInfo),") in sample information file does not match number of sector files (",length(seqFilePaths),") found in the folder.")
    }
    
    return(SimpleList("SequencingFolderPath"=SequencingFolderPath, "seqFilePaths"=seqFilePaths, "seqfilePattern"=seqfilePattern, "sampleInfoFilePath"=sampleInfoFilePath, "sectors"=sampleInfo, "callHistory"=match.call()))
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
#'  	\item sector => region/quadrant of the sequencing plate the sample comes from. If files have been split by samples apriori, then the filename associated per sample without the extension. If this is a filename, then be sure to enable 'alreadyDecoded' parameter in \code{\link{findBarcodes}} or \code{\link{decodeByBarcode}}, since contents of this column is pasted together with 'seqfilePattern' parameter in \code{\link{read.SeqFolder}} to find the appropriate file needed.
#'  	\item barcode => unique 4-12bp DNA sequence which identifies the sample. If providing filename as sector, then leave this blank.
#'  	\item primerltrsequence => DNA sequence of the viral LTR primer with viral LTR sequence following the primer landing site. If already trimmed, then mark this as SKIP.
#'  	\item samplename => Name of the sample associated with the barcode
#'  	\item sampledescription => Detailed description of the sample
#'  	\item gender => sex of the sample: male or female
#'  	\item species => species of the sample: homo sapien, mus musculus, etc.
#'  	\item freeze => UCSC freeze to which the sample should be aligned to.
#'  	\item linkersequence => DNA sequence of linker adaptor to found following the genomic sequence. If already trimmed, then mark this as SKIP.
#'  	\item restrictionenzyme => Restriction enzyme used for digestion and sample recovery. Can also be Fragmentase or Sonication!
#' 		}
#'  \item Metadata Parameter Column Description:
#'   \itemize{
#'  	\item ltrbitsequence => DNA sequence of the viral LTR following the primer landing site.
#'  	\item ltrbitidentity => percent of LTR bit sequence to match during the alignment. Default is 1.
#'  	\item primerltridentity => percent of primer to match during the alignment. Default is .85
#'  	\item linkeridentity => percent of linker sequence to match during the alignment. Default is 0.55. Only applies to non-primerID/random tag based linker search.
#'  	\item primeridinlinker => whether the linker adaptor used has primerID/random tag in it? Default is FALSE.
#'  	\item primeridinlinkeridentity1 => percent of sequence to match before the random tag. Default is 0.75. Only applies to primerID/random tag based linker search and when primeridinlinker is TRUE.
#'  	\item primeridinlinkeridentity2 => percent of sequence to match after the random tag. Default is 0.50. Only applies to primerID/random tag based linker search and when primeridinlinker is TRUE.
#'  	\item celltype => celltype information associated with the sample
#'  	\item user => name of the user who prepared or processed the sample
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
read.sampleInfo <- function(sampleInfoPath=NULL, splitBySector=TRUE, interactive=TRUE) {
    ## read file and make sampleInfo object with sample to metadata associations ##
    if(is.null(sampleInfoPath)) {
        stop("No sample information file path provided.")
    }
    
    sampleInfoPath <- normalizePath(sampleInfoPath,mustWork=TRUE)
    
    requiredCols <- c('sector', 'barcode', 'primerltrsequence', 'samplename', 'sampledescription', 'gender', 'species', 'freeze', 'linkersequence', 'restrictionenzyme')
    metaDataCols <- c('ltrbitsequence'='', 'ltrbitidentity'=1, 'primerltridentity'=.85, 'linkeridentity'=.55, 'primeridinlinker'=FALSE, 'primeridinlinkeridentity1'=.75, 'primeridinlinkeridentity2'=.50, 'celltype'='', 'user'=Sys.getenv("USER"), 'startwithin'=3, 'alignratiothreshold'=.7, 'clustersiteswithin'=5, 'keepmultihits'=TRUE , 'genomicpercentidentity'=.98, 'processingdate'=format(Sys.time(), "%Y-%m-%d "))
    
    if(grepl('.xls.?$',sampleInfoPath)) {
        sampleInfo <- unique(read.xlsx(sampleInfoPath,sheetIndex=1,stringsAsFactors=FALSE))
    } else {
        sampleInfo <- unique(read.delim(sampleInfoPath,stringsAsFactors=FALSE))
    }
    names(sampleInfo) <- tolower(gsub("\\.|-|_","",names(sampleInfo)))
    
    # check for required columns
    ColsNotThere <- !requiredCols %in% names(sampleInfo)
    if (any(ColsNotThere)) {
        absentCols <- requiredCols[ColsNotThere]
        info <- paste("Following required column(s) is absent from the Sample Info file:",paste(absentCols,sep="",collapse=", "),sep="")
        stop(info)
    }
    
    # add missing meta data columns
    metaColsNotThere <- !names(metaDataCols) %in% names(sampleInfo)
    if(any(metaColsNotThere)) {
        sampleInfo <- cbind(sampleInfo,as.data.frame(t(metaDataCols[metaColsNotThere]),stringsAsFactors = FALSE))        
    }
    
    # do some formatting to avoid later hassels!
    for(column in c('sector', 'barcode', 'primerltrsequence', 'ltrbitsequence', 'samplename', 'linkersequence', 'restrictionenzyme')) {
        sampleInfo[,column] <- gsub(" ","",sampleInfo[,column])
        if(column %in% c('barcode', 'primerltrsequence', 'ltrbitsequence', 'linkersequence', 'restrictionenzyme')) {
            sampleInfo[,column] <- toupper(sampleInfo[,column])
        }
    }
    
    # confirm ltr bit is correct
    ltrbitTest <- sampleInfo$primerltrsequence=="SKIP"
    if(any(ltrbitTest)) { ## add SKIP to ltrbit as well if primerltrsequence has been trimmed
    	tofix <- which(ltrbitTest)
    	message("adding SKIP to ltrbitsequence to ",length(tofix)," sample since primerltrsequence has been trimmed.")
    	sampleInfo$ltrbitsequence[tofix] <- "SKIP"
    }
    
    ltrbitTest <- nchar(sampleInfo$ltrbitsequence)==0 | sampleInfo$ltrbitsequence==""
    if(any(ltrbitTest)) {
        tofix <- which(ltrbitTest)
        if(interactive) {
            message("LTR bit not found for ",length(tofix)," samples. Use last 7 bases of the LTR primer as the LTR bit? (y/n)")
            choice <- scan(what=character(0),n=1,quiet=TRUE,multi.line=FALSE)
        } else {
            message("LTR bit not found for ",length(tofix)," samples. Using last 7 bases of the LTR primer as the LTR bit.")
            choice <- "y"
        }
        if(tolower(choice)=="y") {
            sampleInfo$ltrbitsequence <- substr(sampleInfo$primerltrsequence,nchar(sampleInfo$primerltrsequence)-6, nchar(sampleInfo$primerltrsequence))
            sampleInfo$primerltrsequence <- substr(sampleInfo$primerltrsequence,1,nchar(sampleInfo$primerltrsequence)-7)
        } else {
            warning("No LTR bit sequence found for following samples: ",paste(sampleInfo$samplename[tofix],sep="",collapse=", "),immediate.=TRUE)
        }        
    }
    
    # check if samplenames are up to the expectations
    samplenametest <- nchar(sampleInfo$samplename)==0 | sampleInfo$samplename==""
    if(any(samplenametest)) {
        stop("No sample names found in following rows of the sample information file ",sampleInfoPath," : ",paste(which(samplenametest),sep="",collapse=", "))
    }
    
    # check for sectors and their usage
    sectortest <- nchar(sampleInfo$sector)==0 | sampleInfo$sector=="" | is.na(sampleInfo$sector)
    if(any(sectortest)) {
        tofix <- which(sectortest)
        if(interactive) {
            message("Sector information not found for ",length(tofix)," samples. Which sector are they from? (1,2,4,etc)")
            choice <- scan(what=character(0),quiet=TRUE,multi.line=FALSE)
        } else {
            message("Sector information not found for ",length(tofix)," samples. Assuming they are from sector 1.")
            choice <- "1"
        }
        if(length(choice)>0) {
            sampleInfo$sector[tofix] <- unlist(strsplit(choice,","))
        } else {
            stop("No Sector information found for following samples: ",paste(sampleInfo$samplename[tofix],sep="",collapse=", "))
        }
    }
    sampleInfo$sector <- gsub("\\.0$","",as.character(sampleInfo$sector))
    
    sampleSectorTest <- table(paste(sampleInfo$samplename,sampleInfo$sector))
    if(any(sampleSectorTest>1)) {
        stop("Duplicate sample names found on same quadrant in the sample information file ",sampleInfoPath," : ",paste(sampleSectorTest[sampleSectorTest>1],sep="",collapse=", "))
    }
    
    # prepare the sample info object!
    if(splitBySector) {            
        sampleInfo <- SimpleList(split(sampleInfo,sampleInfo$sector))
        for(sector in 1:length(sampleInfo)) { 
            sampleData <- SimpleList(split(sampleInfo[[sector]],sampleInfo[[sector]]$samplename))
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
    stopifnot(class(dnaSet)=="DNAStringSet")
    if(is.null(names(dnaSet))) {
        message("No names attribute found in dnaSet object...using artifically generated names")
        names(dnaSet) <- paste("read",1:length(dnaSet),sep="-")
    }
    dnaSet <- sort(dnaSet)
    ranks <- rank(dnaSet)
    counts <- table(ranks)
    isDuplicate <- duplicated(ranks)
    seqToRank <- data.frame(ranks,counts=as.numeric(counts[as.character(ranks)]),row.names=names(dnaSet),stringsAsFactors = FALSE)
    seqToRank <- seqToRank[!isDuplicate,]
    dnaSet <- dnaSet[!isDuplicate]    
    names(dnaSet) <- paste(names(dnaSet),"counts=",seqToRank[names(dnaSet),"counts"],sep="")
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
replicateReads <- function(dnaSet,counts=NULL) {
    stopifnot(class(dnaSet)=="DNAStringSet")
    if(is.null(counts)) {
        if(is.null(names(dnaSet))) {
            stop("No names attribute found in dnaSet object")
        }
        counts <- as.numeric(sub(".+counts=(\\d+)","\\1",names(dnaSet)))
        if(all(is.na(counts))) {
            stop("No counts=X marker found at the end of definition line or names attribute in dnaSet object")
        }
    }
    if (length(counts)==1) {
        counts <- rep(counts,length(dnaSet))
    }
    
    ids <- unlist(sapply(counts,function(x) 1:x))
    deflines <- paste(rep(sub("(.+)counts=.+","\\1",names(dnaSet)),times=counts),ids,sep="")
    dnaSet <- rep(dnaSet,times=counts)
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
removeReadsWithNs <- function(dnaSet,maxNs=5,consecutive=TRUE) {
    if(consecutive) {
        good.row <- grepl(paste(rep("N",maxNs),collapse=""),dnaSet,fixed=TRUE)
    } else {
        res <- alphabetFrequency(dnaSet)
        good.row <- res[,"N"] <= maxNs
    }
    return(dnaSet[good.row])
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
#' #splitByBarcode(c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2"),dnaSet,showStats=TRUE,returnUnmatched=TRUE)
#'
splitByBarcode <- function(barcodesSample, dnaSet, trimFrom=NULL, showStats=FALSE, returnUnmatched=FALSE, dereplicate=FALSE) {
    if(is.null(barcodesSample) | length(barcodesSample)==0) {
        stop("No barcodes to samples association vector provided in parameter barcodesSample.")
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
    deflines <- sub("^(\\S+) .+$","\\1",names(dnaSet)[good.row],perl = T)
    
    ## if primer bases were utilized for tiebreaker, use the original length instead of modified for trimming.
    if(is.null(trimFrom)) {        
        trimFrom <- barcodelen+1
    }
    
    ## remove sequences with unknown barcode and trim barcode itself ##
    unmatched <- DNAStringSet(dnaSet[!good.row])
    dnaSet <- DNAStringSet(dnaSet[good.row],start=trimFrom)
    names(dnaSet) <- deflines
    
    if(showStats) {
        message("Number of Sequences with no matching barcode: ",as.numeric(table(good.row)['FALSE']))
        message("Number of Sequences decoded:")
        print(as.data.frame(table(sampleNames)))
    }
        
    dnaSet <- as.list(split(dnaSet,as.character(sampleNames)))
    
    if(dereplicate) {
        message("Dereplicating reads.")
        dnaSet <- lapply(dnaSet,dereplicateReads)
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
#' @param sampleInfo sample information SimpleList object created using \code{\link{read.sampleInfo}}, which holds barcodes and sample names per sector/quadrant or a character vector of barcodes to sample name associations. Ex: c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2",...)
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
decodeByBarcode <- function(sampleInfo, sector=NULL, dnaSet=NULL, showStats=FALSE, returnUnmatched=FALSE, dereplicate=FALSE, alreadyDecoded=FALSE) {
    ## tried PDict()...and its slower than this awesome code! ##    
    if(class(sampleInfo)=="SimpleList") {
        if(is.null(sector)) {
            stop("No sector provided in parameter sector.")
        }
        sectors <- sector <- as.character(sector)
        if(length(sectors)==1 & tolower(sectors)=="all") {
            sectors <- names(sampleInfo$sectors)
        }
        if(any(!sectors %in% names(sampleInfo$sectors))) {
            stop("Following sectors not found in names(sampleInfo$sectors): ",sectors[!sectors %in% names(sampleInfo$sectors)])
        }
                
        for(sector in sectors) {
        	## check everything is cool with the provided barcodes first before reading the sequences! ##
            message("Decoding sector: ",sector)            
			## prepare a vector of barcode to sample associations ##
			sampleBarcodes <- toupper(extractFeature(sampleInfo,sector=sector,feature="barcode")[[sector]])
			barcodesSample <- structure(names(sampleBarcodes), names=as.character(sampleBarcodes))
	
			if (length(table(nchar(as.character(sampleBarcodes))))>1) {
				stop("Multiple barcode lengths found.")
			}
			
			## length of barcodes before any modifications done later if any! ##
			realbarcodelen <- unique(nchar(as.character(sampleBarcodes)))
			
			if (any(table(as.character(sampleBarcodes))>1)) {
				message("Duplicate barcode found on this sector.\nPlease choose from one of the options below\n 1: Pick first few bases of primer for tiebreaker? (This could be dangerous if the sequencing run has too many errors!)\n 2: Use the last sample associated with the duplicate as the primary sample?\n 3: Do not do anything.")
				
				choice <- scan(what=integer(0),n=1,quiet=TRUE,multi.line=FALSE)
				
				if(choice==1) {  
					message("Enter # of bases to use from primer:")
					howmany <- scan(what=integer(0),n=1,quiet=TRUE,multi.line=FALSE)
					samplePrimers <- toupper(extractFeature(sampleInfo,sector=sector,feature="primerltrsequence")[[sector]])
					newBarcodes <- toupper(paste(sampleBarcodes,substr(samplePrimers,1,howmany),sep=""))                                                    
					if (any(table(newBarcodes)>1)) {
						stop("Tiebreaking failed...try choosing high number of bases from primer??? \n",paste(names(which(table(newBarcodes)>1)),collapse=", "))
					}
					barcodesSample <- structure(names(sampleBarcodes), names=newBarcodes)
				} else if(choice==2) {
					message("Overwriting duplicate samples associated with the same barcode...")
				} else {
					stop("Aborting due to duplicate barcode found on this sector")
				}
			}
			            
            dnaSet <- read.seqsFromSector(sampleInfo,sector)
            
            if(alreadyDecoded) {
            	if(length(barcodesSample)>1) {
            		stop("alreadyDecoded parameter is set to TRUE. There shouldn't be more than one sample associated to a sequence file.")
            	}
				names(dnaSet) <- sub("^(\\S+) .+$","\\1",names(dnaSet),perl = T)
				dnaSet <- as.list(split(dnaSet,rep(as.character(barcodesSample),length(dnaSet))))
            } else {
            	dnaSet <- splitByBarcode(barcodesSample, dnaSet, trimFrom=realbarcodelen+1, showStats=showStats, returnUnmatched=returnUnmatched, dereplicate=dereplicate)
            }
            
            for(samplename in names(dnaSet)) {
                if(samplename=="unDecodedSeqs") {                                        
                    metadata(sampleInfo$sectors[[sector]]) <- append(metadata(sampleInfo$sectors[[sector]]), list("unDecodedSeqs"=dnaSet[[samplename]]))            
                } else {
                    sampleInfo$sectors[[sector]]$samples[[samplename]]$decoded <- dnaSet[[samplename]]
                }
            }
            metadata(sampleInfo$sectors[[sector]]) <- append(metadata(sampleInfo$sectors[[sector]]), list("decodedBy"=barcodesSample))
        }
        sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
        decoded <- sampleInfo
        cleanit <- gc()
    } else {
        decoded <- splitByBarcode(sampleInfo, dnaSet, trimFrom=NULL, showStats=showStats, returnUnmatched=returnUnmatched, dereplicate=dereplicate)
        cleanit <- gc()
    }
    return(decoded)
}

#' @export
findBarcodes <- decodeByBarcode

#' Align a short pattern to variable length target sequences.
#'
#' Align a fixed length short pattern sequence (i.e. primers or adaptors) to subject sequences using \code{\link{pairwiseAlignment}}. This function uses default of type="overlap", gapOpening=-1, and gapExtension=-1 to align the patternSeq against subjectSeqs. One can adjust these parameters if prefered, but not recommended. This function is meant for aligning a short pattern onto large collection of subjects. If you are looking to align a vector sequence to subjects, then please use BLAT.
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
#' @param ... extra parameters for \code{\link{pairwiseAlignment}}
#'
#' @note For qualityThreshold, the alignment score is calculated by (matches*2)-(mismatches+gaps) which programatically translates to round(nchar(patternSeq)*qualityThreshold)*2. Gaps and mismatches are weighed equally with value of -1. If qualityThreshold is 1, then its is a full match, if 0, then any match is accepted which is useful in searching linker sequences at 3' end. Beware, this function only searches for the pattern sequence in one orientation. If you are expecting to find the pattern in both orientation, you might be better off using BLAST!
#'
#' @return IRanges object with starts, stops, and names of the aligned sequences. If returnLowScored or returnUnmatched = T, then a CompressedIRangesList where x[["hits"]] has the good scoring hits, x[["Rejected"]] has the failed to match qualityThreshold hits, and x[["Absent"]] has the hits where the aligned bit is <=10\% match to the patternSeq.
#'
#' @seealso \code{\link{primerIDAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{doRCtest}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #pairwiseAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE)
#' #pairwiseAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE,qualityThreshold=0.5)
#'
pairwiseAlignSeqs <- function(subjectSeqs=NULL, patternSeq=NULL, side="left", qualityThreshold=1, showStats=FALSE, bufferBases=5, doRC=TRUE, returnUnmatched=FALSE, returnLowScored=FALSE, ...) {
    if(is.null(subjectSeqs) | is.null(patternSeq)) {
        stop("subjectSeqs/patternSeq is empty. Please supply reads to be aligned")
    }
    
    if(length(patternSeq)>1) {
        stop("More than 1 patternSeq is defined. Please only supply one pattern.")
    }
    
    ## give names if not there for troubleshooting purpose in later steps
    if(is.null(names(subjectSeqs))) {
        names(subjectSeqs) <- paste("read",1:length(subjectSeqs))
    }
    
    ## only get the relevant side of subject sequence with extra bufferBases to account for indels/mismatches & save memory while searching & avoid searching elsewhere in the sequence
    if(tolower(side)=="left") {
        badSeqs <- DNAStringSet()
        culprits <- width(subjectSeqs) < (nchar(patternSeq)+bufferBases)
        if(any(culprits)) {
            badSeqs <- subjectSeqs[culprits]
            message(length(badSeqs)," sequences were removed from aligning since they were shorter than pattern getting aligned: ",(nchar(patternSeq)+bufferBases),"bp")            
            subjectSeqs <- subjectSeqs[!culprits]            
        }
        subjectSeqs2 <- subseq(subjectSeqs,start=1,end=(nchar(patternSeq)+bufferBases))
        overFromLeft <- rep(0,length(subjectSeqs))
    } else if (tolower(side)=="right") { 
        overFromLeft <- width(subjectSeqs)-(nchar(patternSeq)+bufferBases)
        overFromLeft[overFromLeft<1] <- 1
        subjectSeqs2 <- subseq(subjectSeqs,start=overFromLeft)
    } else {
        subjectSeqs2 <- subjectSeqs
        overFromLeft <- rep(0,length(subjectSeqs))
    }
    
    ## search both ways to test which side yields more hits!
    if(doRC) {
        patternSeq <- doRCtest(subjectSeqs2, patternSeq, qualityThreshold)            
    }
    
    ## do not split this into a multicore function since its faster to align all sequences at once then breaking into smaller chunks
    ## type=overlap is best for primer trimming...see Biostrings Alignment vignette
    if(any(names(match.call()) %in% c("type","gapOpening","gapExtension"))) {
    	hits <- pairwiseAlignment(subjectSeqs2, patternSeq, ...)        
    } else {
    	hits <- pairwiseAlignment(subjectSeqs2, patternSeq, type="overlap", gapOpening=-1, gapExtension=-1, ...)
    }
    
    stopifnot(length(hits)==length(subjectSeqs2))
    
    scores <- round(score(hits))
    highscored <- scores >= round(nchar(patternSeq)*qualityThreshold)*2
    unmatched <- nchar(hits) <= round(nchar(patternSeq)*.1) ## basically a small subset of highscored
        
    # no point in showing stats if all sequences are a potential match! #
    if(showStats & qualityThreshold!=0) {
        message("Total of ",as.numeric(table(highscored)['FALSE'])," did not have the defined pattern sequence (",patternSeq,") that passed qualityThreshold on the ",side," side")
    }

    ## extract starts-stops of the entire pattern hit ##
    starts <- start(pattern(hits))
    ends <- end(pattern(hits))
    namesq <- names(subjectSeqs)
    hits <- IRanges(start=starts+overFromLeft-ifelse(side=="right",2,0), end=ends+overFromLeft-ifelse(side=="right",2,0), names=namesq)
    rm("scores","subjectSeqs2","subjectSeqs","starts","ends","namesq")
    
    ## no need to test if there were any multiple hits since pairwiseAlignment will only output one optimal alignment...see the man page.
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
        hits <- hitstoreturn; rm(hitstoreturn)
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
#' @note For qualityThreshold1 & qualityThreshold2, the alignment score is calculated by (matches*2)-(mismatches+gaps) which programatically translates to round(nchar(patternSeq)*qualityThreshold)*2. Gaps and mismatches are weighed equally with value of -1. If qualityThreshold is 1, then its is a full match, if 0, then any match is accepted which is useful in searching linkers at 3' end. Beware, this function only searches for the pattern sequence in one orientation. If you are expecting to find the pattern in both orientation, you might be better off using BLAST!
#'
#' @return A CompressedIRangesList of length two, where x[["hits"]] is hits covering the entire patternSeq, and x[["primerIDs"]] is the potential primerID region. If returnUnmatched = T, then x[["Absent"]] is appended which includes reads not matching the first part of patternSeq. If returnRejected=TRUE, then x[["Rejected"]] includes reads that only matched one part of patternSeq or places where no primerID was found in between two part of patternSeq, and x[["RejectedprimerIDs"]] includes primerIDs that didn't match the correct length. If doAnchored=TRUE, then x[["unAnchoredprimerIDs"]] includes reads that didn't match the base before and after primer ID on patternSeq.
#'
#' @seealso \code{\link{vpairwiseAlignSeqs}}, \code{\link{pairwiseAlignSeqs}}, \code{\link{doRCtest}}
#'
#' @export
#'
#' @examples 
#' #primerIDAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE)
#' #primerIDAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE,qualityThreshold1=0.5)
#'
primerIDAlignSeqs <- function(subjectSeqs=NULL, patternSeq=NULL, qualityThreshold1=0.75, qualityThreshold2=0.50, doAnchored=FALSE, doRC=TRUE, returnUnmatched=FALSE, returnRejected=FALSE, showStats=FALSE, ...) {
    if(is.null(subjectSeqs) | is.null(patternSeq)) {
        stop("subjectSeqs/patternSeq is empty. Please supply reads to be aligned")
    }
    
    if(length(patternSeq)>1) {
        stop("More than 1 patternSeq is defined. Please only supply one pattern.")
    }
    
    ## give names if not there for troubleshooting purpose in later steps
    if(is.null(names(subjectSeqs))) {
        names(subjectSeqs) <- paste("read",1:length(subjectSeqs))
    }

    ## make sure there are Ns in the patternSeq for considering primerIDs
    if(length(unlist(gregexpr("N",patternSeq)))<4) {
        stop("There should be minimum of atleast 4 Ns in patternSeq to be considered as a primerID sequence.")
    }
    
    ## get the right orientation of the supplied patternSeq to peform proper search at later two step search phase. 
    if(doRC) {
        patternSeq <- doRCtest(subjectSeqs, patternSeq)            
    }
    
    primerIDpos <- unlist(gregexpr("N",patternSeq))
    
    ## perform primerID extraction by breaking the pattern into two parts for sanity sakes due to homopolymers ##
    pattern1 <- as.character(subseq(DNAString(patternSeq),1,primerIDpos[1]-1))
    pattern2 <- as.character(subseq(DNAString(patternSeq),primerIDpos[length(primerIDpos)]+1))

    pattern1.hits <- pairwiseAlignSeqs(subjectSeqs, pattern1, "middle", qualityThreshold=qualityThreshold1, doRC=FALSE, returnUnmatched=TRUE, ...)
    pattern2.hits <- pairwiseAlignSeqs(subjectSeqs, pattern2, "middle", qualityThreshold=qualityThreshold2, doRC=FALSE, ...)
    
    ## set aside reads which has no match to the pattern1...most likely mispriming if primerID is on 5' or read was too loong if on 3'
    ## hits returned from pairwiseAlignSeqs will be filtered for low scored hits...so no need to check for those from subjectSeqs
    unmatched <- pattern1.hits[["Absent"]]
    pattern1.hits <- pattern1.hits[["hits"]]
    
    ## remove reads which only have a match to one of the patterns...crossover most likely!
    rejected1 <- setdiff(names(pattern1.hits),names(pattern2.hits))
    rejected2 <- setdiff(names(pattern2.hits),names(pattern1.hits))
    rejected <- c(pattern1.hits[names(pattern1.hits) %in% rejected1], pattern2.hits[names(pattern2.hits) %in% rejected2])
    rm("rejected1","rejected2")
    if(showStats) { message("Removed ",length(rejected)," read(s) for only matching one of pattern1 or pattern2") }
    
    ## use only reads which match to both sides of the patterns.
    good.rows <- intersect(names(pattern1.hits),names(pattern2.hits))
    pattern1.hits <- pattern1.hits[names(pattern1.hits) %in% good.rows]
    pattern2.hits <- pattern2.hits[names(pattern2.hits) %in% good.rows]
    
    stopifnot(identical(names(pattern1.hits),names(pattern2.hits)))
    stopifnot(identical(names(pattern1.hits),good.rows))
    ## make sure there is no overlap of ranges between pattern1.hits & pattern2.hits
    ## if there is, then no primerID was found...remove it
    badAss <- end(pattern1.hits) >= start(pattern2.hits)
    if(any(badAss)) {
        rejected <- c(rejected, pattern1.hits[badAss], pattern2.hits[badAss])
        pattern1.hits <- pattern1.hits[!badAss]
        pattern2.hits <- pattern2.hits[!badAss]
        good.rows <- good.rows[!badAss]
        message("Removed ",table(badAss)["TRUE"]," read(s) for not having primerID present between pattern1 and pattern2")
    }
    
    hits <- IRanges(start=start(pattern1.hits), end=end(pattern2.hits), names=good.rows)        
    primerIDs <- IRanges(start=end(pattern1.hits)+1, end=start(pattern2.hits)-1, names=good.rows)        
        
    if(length(hits)==0) {
        stop("No hits found that matched both sides of patternSeq with primerID in the middle.")
    }    
    
    ## do anchored search for only sequences that matched both sides of patternSeq
    if(doAnchored) {
        message("Found ",length(primerIDs)," total primerIDs before anchored filter.")

        ## get anchors of bases flanking Ns
        anchorBase.s <- substr(patternSeq,primerIDpos[1]-1,primerIDpos[1]-1)
        anchorBase.e <- substr(patternSeq,primerIDpos[length(primerIDpos)]+1,primerIDpos[length(primerIDpos)]+1)
        
        anchorBase.s.i <- trimSeqs(subjectSeqs,resize(pattern1.hits,width=1,fix="end"))
        anchorBase.e.i <- trimSeqs(subjectSeqs,resize(pattern2.hits,width=1,fix="start"))
        rows <- anchorBase.s==as.character(anchorBase.s.i) & anchorBase.e==as.character(anchorBase.e.i)
        
        unAnchored <- hits[!rows]
        primerIDs <- primerIDs[rows]
        hits <- hits[rows]
        message("Found ",length(primerIDs)," total primerIDs after anchored filter.")
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
        message("Removed ",table(badAss)["TRUE"]," read(s) for not having right primerID length")
    }

    hits <- IRangesList("hits"=hits,"primerIDs"=primerIDs)
    
    if(exists("unAnchored")) {
        if(length(unAnchored)>0) { hits <- append(hits,IRangesList("unAnchoredprimerIDs"=unAnchored)) }
    }
    
    if(returnUnmatched & length(unmatched)>0) {
        hits <- append(hits,IRangesList("Absent"=unmatched))
    }
    
    if(returnRejected) {
        if(length(rejected)>0) { hits <- append(hits,IRangesList("Rejected"=rejected)) }
        if(exists("rejectedprimerIDs")) { 
            if(length(rejectedprimerIDs)>0) { hits <- append(hits,IRangesList("RejectedprimerIDs"=rejectedprimerIDs)) } 
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
#' @param ... extra parameters for \code{\link{vmatchPattern}} except for 'max.mismatch' since it's calculated internally using the 'qualityThreshold' parameter.
#'
#' @note Beware, this function only searches for the pattern sequence in one orientation. If you are expecting to find the pattern in both orientation, you might be better off using BLAST!
#'
#' @return IRanges object with starts, stops, and names of the aligned sequences.
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{doRCtest}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #vpairwiseAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE)
#' #vpairwiseAlignSeqs(subjectSeqs,patternSeq,showStats=TRUE,qualityThreshold=0.5)
#'
vpairwiseAlignSeqs <- function(subjectSeqs=NULL, patternSeq=NULL, side="left", qualityThreshold=1, showStats=FALSE, bufferBases=5, doRC=TRUE, ...) {
    if(is.null(subjectSeqs) | is.null(patternSeq)) {
        stop("subjectSeqs/patternSeq is empty. Please supply reads to be aligned")
    }

    ## give names if not there for troubleshooting purpose in later steps
    if(is.null(names(subjectSeqs))) {
        names(subjectSeqs) <- paste("read",1:length(subjectSeqs))
    }
    
    ## only get the relevant side of subject sequence with extra bufferBases to account for indels/mismatches & save memory while searching & avoid searching elsewhere in the sequence
    if(tolower(side)=="left") {
        badSeqs <- DNAStringSet()
        culprits <- width(subjectSeqs) < (nchar(patternSeq)+bufferBases)
        if(any(culprits)) {
            badSeqs <- subjectSeqs[culprits]
            message(length(badSeqs)," sequences were removed from aligning since they were shorter than pattern getting aligned: ",(nchar(patternSeq)+bufferBases),"bp")
            subjectSeqs <- subjectSeqs[!culprits]            
        }
        subjectSeqs2 <- subseq(subjectSeqs,start=1,end=(nchar(patternSeq)+bufferBases))
        overFromLeft <- rep(0,length(subjectSeqs))
    } else if (tolower(side)=="right") { 
        overFromLeft <- width(subjectSeqs)-(nchar(patternSeq)+bufferBases)
        overFromLeft[overFromLeft<1] <- 1
        subjectSeqs2 <- subseq(subjectSeqs,start=overFromLeft)
    } else {
        subjectSeqs2 <- subjectSeqs
        overFromLeft <- rep(0,length(subjectSeqs))
    }
    
    ## search both ways to test which side yields more hits!        
    if(doRC) {        
        patternSeq <- doRCtest(subjectSeqs2, patternSeq, qualityThreshold)            
    }
    
    ## do not split, vmatchPattern works faster when given one contiguous dnastringset
    hits <- vmatchPattern(patternSeq, subjectSeqs2, max.mismatch=round(nchar(patternSeq)*(1-qualityThreshold)), ...)
    hits <- unlist(hits, recursive=TRUE, use.names=TRUE)
        
    ## test if there were any multiple hits which are overlapping and if they are reduce them
    counts <- Rle(names(hits))
    if(any(runLength(counts)>1)) {
        reduced <- reduce(GRanges(seqnames=names(hits),IRanges(start=start(hits),end=end(hits))))
        counts <- seqnames(reduced)
        hits <- ranges(reduced)
        names(hits) <- as.character(seqnames(reduced))
        rm(reduced)
        if(any(runLength(counts)>1)) {
            message(paste("More than 1 pattern (",patternSeq,") match found for:",paste(runValue(counts)[runLength(counts)>1],collapse=",")))
            message("\nUsing the latter occuring hit as the dominant for each read.")
            toremove <- c()
            for(culprits in as.character(runValue(counts)[runLength(counts)>1])) {
                rows <- which(names(hits) %in% culprits)
                toremove <- c(toremove,rows[1:length(rows)-1])
            }
            hits <- hits[-toremove]
            counts <- Rle(names(hits))
            if(any(runLength(counts)>1)) {
                stop(paste("More than 1 pattern unresolved (",patternSeq,") match found for:",paste(runValue(counts)[runLength(counts)>1],collapse=",")))
            }
        }
    }
        
    good.row <- names(subjectSeqs2) %in% names(hits)
    
    if(showStats) {
        message("Total of ",as.numeric(table(good.row)['FALSE'])," did not have the defined pattern sequence (",patternSeq,") that passed qualityThreshold on the ",side," side")
    }  
    
    starts <- start(hits)
    ends <- end(hits)
    namesq <- names(hits)
    rm("hits","subjectSeqs2")     
    
    hits <- IRanges(start=starts+overFromLeft[good.row]-ifelse(side=="right",2,0), end=ends+overFromLeft[good.row]-ifelse(side=="right",2,0), names=namesq)
    
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
doRCtest <- function(subjectSeqs=NULL, patternSeq=NULL, qualityThreshold=0.5, core.use=2) {
    if(is.null(subjectSeqs) | is.null(patternSeq)) {
        stop("subjectSeqs/patternSeq is empty. Please supply reads to be aligned")
    }

    if(core.use==1) { registerDoSEQ() }
    
    patternSeq.rc <- as.character(reverseComplement(DNAString(patternSeq)))
    hits <- foreach(x=iter(c(patternSeq,patternSeq.rc)), .inorder=TRUE, .errorhandling="pass", .export=c("subjectSeqs","qualityThreshold"), .packages="Biostrings") %dopar% {
        counts <- pmin(vcountPattern(x, subjectSeqs, max.mismatch=round(nchar(x)*(1-qualityThreshold))),1)
        sum(counts)
    }
    
    if(all(hits==0)) {
        stop("No hits found")
    }
        
    if(hits[[1]] < hits[[2]]) {
        message("There were less/no good hits found for pattern ",patternSeq,"\nthan its reverse complement...",patternSeq.rc,"\nUsing the latter to perform the searching.")
        patternSeq <- patternSeq.rc
    }
    
    return(patternSeq)
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
#' @param ... extra parameters to be passed to either \code{\link{vpatternMatch}} or \code{\link{pairwiseAlignment}} depending on 'alignWay' parameter.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new data added under each sector and sample. New data attributes include: primed
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{extractFeature}}, \code{\link{extractSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}}, \code{\link{findLinkers}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #findPrimers(sampleInfo,showStats=TRUE)
#' #findPrimers(sampleInfo,alignWay="slow",showStats=TRUE)
#'
findPrimers <- function(sampleInfo, alignWay="slow", showStats=FALSE, doRC=FALSE, parallel=TRUE, samplenames=NULL, ...) {    
    stopifnot(class(sampleInfo)=="SimpleList")
    
    if(!parallel) { registerDoSEQ() }
    
    ## test if there are decoded sequences in the sampleinfo object ##
    decoded <- extractFeature(sampleInfo,feature="decoded")
    samplesDecoded <- sapply(decoded,names,simplify=FALSE)
    sectorsDecoded <- names(which(sapply(sapply(decoded,length),">",0)))
    rm(decoded)
    gc()

    if(length(sectorsDecoded)==0) {
        stop("No decoded information found in sampleInfo...did you run decodeByBarcode()?")
    }
    
    for(sector in sectorsDecoded) {
        message("Processing sector ",sector)
        
        ## prepare sample to primer association ##
        ltrPrimers <- toupper(extractFeature(sampleInfo,sector=sector,feature="primerltrsequence")[[sector]])
        skippers <- ltrPrimers=="SKIP"
        if(!all(skippers) & (length(ltrPrimers[!skippers])==0 | mean(nchar(ltrPrimers[!skippers]))<=5 | all(is.na(ltrPrimers[!skippers])))) {
            stop("Either the primer size is too short (<=5) or no primers are found in sample information object.")
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
            sampleInfo <- addFeature(sampleInfo, sector, skip.samples, 
            						feature="primed", value=structure(rep("SKIPPED",length(skip.samples)), names=skip.samples))
        }
        
        ## dont bother searching if no samples are to be processed! ##
        if(length(samplesToProcess)>0) {			
        	## get the decoded reads ##
			decoded <- extractSeqs(sampleInfo,sector,samplesToProcess,feature="decoded")[[sector]]
		
			## trim the primers ##
			message("\tFinding Primers.")
			primerIdentity <- extractFeature(sampleInfo,sector=sector,feature="primerltridentity")[[sector]]
			stopifnot(length(primerIdentity)>0)
			
			primerTrimmed <- foreach(x=iter(samplesToProcess), .inorder=TRUE, .errorhandling="pass", .export=c("ltrPrimers","primerIdentity","doRC","alignWay", "vpairwiseAlignSeqs", "pairwiseAlignSeqs","decoded"), .packages="Biostrings") %dopar% {
				switch(alignWay,
					fast = vpairwiseAlignSeqs(decoded[[x]], ltrPrimers[[x]], "left", qualityThreshold=(primerIdentity[[x]]-.05), doRC=doRC, ...),
					slow = pairwiseAlignSeqs(decoded[[x]], ltrPrimers[[x]], "left", qualityThreshold=(primerIdentity[[x]]), doRC=doRC, ...)                
				)        
			}
			names(primerTrimmed) <- samplesToProcess
			
			## check if any error occured during alignments ##
			if(any(grepl("simpleError",primerTrimmed))) {
				stop("Error encountered in LTR Trimming function",paste(names(primerTrimmed[grepl("simpleError",primerTrimmed)]),collapse=", "))
			}        
				
			## remove samples with no primer hits from further processing ##
			culprits <- grep("No hits found",primerTrimmed)
			if(length(culprits)>0) {
				message("Following sample(s) had no hits for primer alignment: ",paste(samplesToProcess[culprits],collapse=", "))
				samplesToProcess <- samplesToProcess[-c(culprits)]
				primerTrimmed <- primerTrimmed[-c(culprits)]
			}
			
			cleanit <- gc()
			
			toprint <- as.data.frame(sapply(primerTrimmed,length)); names(toprint) <- "Total"            
			counts <- sapply(decoded,length)
			toprint$PercOfDecoded <- 100*(toprint$Total/counts[rownames(toprint)])
			toprint$SampleName <- rownames(toprint)
			rownames(toprint) <- NULL
	
			if (showStats) {            
				message("Sequence Stats after primer alignment:")
				print(toprint)        
			}
	
			## if <= 5% of sequences found primers...then something is wrong with the primer sequences provided???
			if(mean(toprint$PercOfDecoded)<=5) {
				stop("Something seems to be wrong with the primers provided for each sample. On average <= 5% of sequences found primer match for the entire run!!!")
			}
			
			## modify metadata attribute, write primer coordinates back to sampleInfo object & trim
			message("Adding primer info back to the object")
			sampleInfo <- addFeature(sampleInfo,sector,names(primerTrimmed),feature="primed",value=primerTrimmed)
			rm(primerTrimmed)
			cleanit <- gc()
		}
    }
    sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
    return(sampleInfo)
}

#' Find the 5' LTRs and add results to SampleInfo object. 
#'
#' Given a sampleInfo object, the function finds 5' LTR following the primer for each sample per sector and adds the results back to the object. This is a specialized function which depends on many other functions shown in 'see also section' to perform specialized trimming of 5' viral LTRs found in the sampleInfo object. The sequence itself is never trimmed but rather coordinates of LTR portion is added to primer coordinates and recorded back to the object and used subsequently by \code{\link{extractSeqs}} function to perform the trimming. This function heavily relies on \code{\link{pairwiseAlignSeqs}}.
#'
#' @param sampleInfo sample information SimpleList object outputted from \code{\link{findPrimers}}, which holds decoded and primed sequences for samples per sector/quadrant along with information of sample to LTR associations.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param doRC perform reverse complement search of the defined pattern/LTR sequence. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param samplenames a vector of samplenames to process. Default is NULL, which processes all samples from sampleInfo object.
#' @param ... extra parameters to be passed to \code{\link{pairwiseAlignment}}.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new data added under each sector and sample. New data attributes include: LTRed
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{extractFeature}}, \code{\link{extractSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findPrimers}}, \code{\link{findLinkers}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #findLTRs(sampleInfo,showStats=TRUE)
#'
findLTRs <- function(sampleInfo, showStats=FALSE, doRC=FALSE, parallel=TRUE, samplenames=NULL, ...) {    
    stopifnot(class(sampleInfo)=="SimpleList")
    
    if(!parallel) { registerDoSEQ() }
    
    ## test if there are primed sequences in the sampleinfo object ##   
    primed <- extractFeature(sampleInfo,feature="primed")
    samplesprimed <- sapply(primed,names,simplify=FALSE)
    sectorsPrimed <- names(which(sapply(sapply(primed,length),">",0)))
    rm(primed)
    gc()

    if(length(sectorsPrimed)==0) {
        stop("No primed information found in sampleInfo...did you run findPrimers()?")
    }
    
    for(sector in sectorsPrimed) {
        message("Processing sector ",sector)
        
        ## prepare sample to LTR bit associations ##
        sampleLTRbits <- toupper(extractFeature(sampleInfo,sector=sector,feature="ltrbitsequence")[[sector]])
        skippers <- sampleLTRbits=="SKIP"
        if(!all(skippers) & (length(sampleLTRbits[!skippers])==0 | mean(nchar(sampleLTRbits[!skippers]))<=1 | all(is.na(sampleLTRbits[!skippers])))) {
            stop("Either LTR bit sequence is too short (<=1) or no LTR bits found in sample information object.")
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
            						feature="LTRed", value=structure(rep("SKIPPED",length(skip.samples)), names=skip.samples))
        }
		
		## dont bother searching if no samples are to be processed! ##
        if(length(samplesToProcess)>0) {
			## get the primer trimmed reads ##
			primerTrimmed <- extractSeqs(sampleInfo,sector,samplesToProcess,feature="primed")[[sector]]
			
			## trim LTRbits using slow method since its the best! ##
			message("\tFinding LTR bits.")
			ltrBitIdentity <- extractFeature(sampleInfo,sector=sector,feature="ltrbitidentity")[[sector]]
			ltrTrimmed <- foreach(x=iter(samplesToProcess), .inorder=TRUE, .errorhandling="pass", .export=c("primerTrimmed","sampleLTRbits","ltrBitIdentity","doRC", "pairwiseAlignSeqs"), .packages="Biostrings") %dopar% {            
				pairwiseAlignSeqs(primerTrimmed[[x]], sampleLTRbits[[x]], "left", qualityThreshold=ltrBitIdentity[[x]], doRC=doRC, ...)
			}
			names(ltrTrimmed) <- samplesToProcess
			
			## check if any error occured during alignments ##
			if(any(grepl("simpleError",ltrTrimmed))) {
				stop("Error encountered in LTR Trimming function",paste(names(ltrTrimmed[grepl("simpleError",ltrTrimmed)]),collapse=", "))
			}
			
			## remove samples with no LTR hits from further processing ##
			culprits <- grep("No hits found",ltrTrimmed)
			if(length(culprits)>0) {
				message("Following sample(s) had no hits for LTR bit alignment: ",paste(samplesToProcess[culprits],collapse=", "))
				samplesToProcess <- samplesToProcess[-c(culprits)]
				ltrTrimmed <- ltrTrimmed[-c(culprits)]
			}    
			
			cleanit <- gc()
	
			toprint <- as.data.frame(sapply(ltrTrimmed,length)); names(toprint) <- "Total"
			counts <- sapply(primerTrimmed,length)
			toprint$PercOfPrimed <- 100*(toprint$Total/counts[rownames(toprint)])        
			toprint$SampleName <- rownames(toprint)
			rownames(toprint) <- NULL    
			
			if (showStats) {                            
				message("Sequence Stats after LTR bit alignment:")
				print(toprint)        
			}
	
			## if <= 5% of sequences found LTRs...then something is wrong with the LTR sequences provided???
			if(mean(toprint$PercOfPrimed)<=5) {
				stop("Something seems to be wrong with the LTRs provided for each sample. On average <= 5% of sequences found LTR match for the entire run!!!")
			}
			
			## modify metadata attribute, add LTR bit coordinates to primer and write back to sampleInfo object & trim
			message("Adding LTR info back to the object")
			for(x in names(ltrTrimmed)) {
				cat('.')
				if(length(ltrTrimmed[[x]])>0) {
					primed.end <- end(sampleInfo$sectors[[sector]]$samples[[x]]$primed[names(ltrTrimmed[[x]])])
					end(ltrTrimmed[[x]]) <- end(ltrTrimmed[[x]]) + primed.end
					start(ltrTrimmed[[x]]) <- start(ltrTrimmed[[x]]) + primed.end
					sampleInfo$sectors[[sector]]$samples[[x]]$LTRed <- ltrTrimmed[[x]]
					rm(primed.end)
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
#' @param ... extra parameters to be passed to \code{\link{pairwiseAlignment}}.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new data added under each sector and sample. New data attributes include: linkered. If linkers have primerID then, primerIDs attribute is appended as well. 
#'
#' @note If no linker matches are found with default options, then try doRC=TRUE. If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}}, \code{\link{findPrimers}}, \code{\link{extractFeature}}, \code{\link{extractSeqs}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #findLinkers(sampleInfo,showStats=TRUE)
#'
findLinkers <- function(sampleInfo, showStats=FALSE, doRC=FALSE, parallel=TRUE, samplenames=NULL, ...) {    
    stopifnot(class(sampleInfo)=="SimpleList")
    
    if(!parallel) { registerDoSEQ() }

    ## test if there are primed or LTRed sequences in the sampleinfo object ##
    primerFlag <- FALSE
    toProcess <- extractFeature(sampleInfo,feature="LTRed")
    toProcessSamples <- sapply(toProcess,names,simplify=FALSE)
    sectors <- names(which(sapply(sapply(toProcess,length),">",0)))
    if(length(sectors)==0) {
        message("No LTRed information found in sampleInfo. Using primed information to find Linkers")
        toProcess <- extractFeature(sampleInfo,feature="primed")
        toProcessSamples <- sapply(toProcess,names,simplify=FALSE)
        sectors <- names(which(sapply(sapply(toProcess,length),">",0)))
        if(length(sectors)==0) {
            stop("No primed information found in sampleInfo...did you run findPrimers()?")
        }
        primerFlag <- TRUE
    }
    rm(toProcess)
    gc()
    
    for(sector in sectors) {
        message("Processing sector ",sector)
        
        ## prepare sample to linker association ##
        sampleLinkers <- toupper(extractFeature(sampleInfo,sector=sector,feature="linkersequence")[[sector]])
        skippers <- sampleLinkers=="SKIP"
        if(!all(skippers) & (length(sampleLinkers[!skippers])==0 | mean(nchar(sampleLinkers[!skippers]))<=10 | all(is.na(sampleLinkers[!skippers])))) {
            stop("Either Linker sequence is too short (<=10) or no Linkers found in sample information object.")
        }
        
        ## get the linker quality thresholds for non primerID based samples ##
        linkerIdentity <- extractFeature(sampleInfo,sector=sector,feature="linkeridentity")[[sector]]
        stopifnot(length(linkerIdentity)>0)

        ## check if any are primerIDed and get their identity thresholds ##
        primerIded <- extractFeature(sampleInfo,sector=sector,feature="primeridinlinker")[[sector]]        
        primerIded.threshold1 <- extractFeature(sampleInfo,sector=sector,feature="primeridinlinkeridentity1")[[sector]]        
        primerIded.threshold2 <- extractFeature(sampleInfo,sector=sector,feature="primeridinlinkeridentity2")[[sector]]        
        
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
            sampleInfo <- addFeature(sampleInfo, sector, skip.samples, 
            						feature="linkered", value=structure(rep("SKIPPED",length(skip.samples)), names=skip.samples))
        }
        
        ## dont bother searching if no samples are to be processed! ##
        if(length(samplesToProcess)>0) {
					
			toProcess <- extractSeqs(sampleInfo,sector,samplesToProcess,feature=ifelse(primerFlag,"primed","LTRed"))[[sector]]
			
			## trim the Linkers ##
			message("\tFinding Linkers.")
			linkerTrimmed <- foreach(x=iter(samplesToProcess), .inorder=TRUE, .errorhandling="pass", .export=c("primerIded","toProcess","sampleLinkers","doRC","linkerIdentity","primerIded.threshold1","primerIded.threshold2","pairwiseAlignSeqs","primerIDAlignSeqs"), .packages="Biostrings") %dopar% {
				if(primerIded[[x]]) {
					primerIDAlignSeqs(toProcess[[x]], sampleLinkers[[x]], doAnchored=TRUE, returnUnmatched=TRUE, returnRejected=TRUE, doRC=doRC, qualityThreshold1=primerIded.threshold1[[x]], qualityThreshold2=primerIded.threshold2[[x]], ...)
				} else {
					pairwiseAlignSeqs(toProcess[[x]], sampleLinkers[[x]], "middle", qualityThreshold=linkerIdentity[[x]], returnUnmatched=TRUE, returnLowScored=TRUE, doRC=doRC, ...) ## use side="middle" since more junk sequence can be present after linker which would fail pairwiseAlignSeqs if side='right'
				}        
			}
			names(linkerTrimmed) <- samplesToProcess
			
			## check if any error occured during alignments ##
			if(any(grepl("simpleError",linkerTrimmed))) {
				stop("Error encountered in Linker Trimming functions",paste(names(linkerTrimmed[grepl("simpleError",linkerTrimmed)]),collapse=", "))
			}
			
			## remove samples with no linker hits from further processing ##
			culprits <- grep("No hits found",linkerTrimmed)
			if(length(culprits)>0) {
				message("Following sample(s) had no hits for Linker alignment: ",paste(samplesToProcess[culprits],collapse=", "))
				samplesToProcess <- samplesToProcess[-c(culprits)]
				linkerTrimmed <- linkerTrimmed[-c(culprits)]
			}        
			
			cleanit <- gc()
			
			toprint <- data.frame("Total"=sapply(sapply(linkerTrimmed,"[[","hits"),length))
			counts <- sapply(toProcess,length)
			toprint$PercOfPrimedOrLTRed <- 100*(toprint$Total/counts[rownames(toprint)]) 
			toprint$SampleName <- rownames(toprint)
			rownames(toprint) <- NULL    
			
			if (showStats) {                            
				message("Sequence Stats after Linker alignment:")
				print(toprint)        
			}
	
			## if <= 5% of sequences found linkers...then something is wrong with the linker sequences provided???
			if(mean(toprint$PercOfPrimedOrLTRed)<=5) {
				stop("Something seems to be wrong with the linkers provided for each sample. On average <= 5% of sequences found linker match for the entire run!!!")
			}
			
			message("Adding linker info back to the object")
			## modify metadata attribute, add linker coordinates to LTRed and write back to sampleInfo object
			## for primerID based samples...write back all the returned attibutes
			for(x in names(linkerTrimmed)) {
				for(y in names(linkerTrimmed[[x]])) {
					cat('.')
					if(primerFlag){
						LTRed.ends <- end(sampleInfo$sectors[[sector]]$samples[[x]]$primed[names(linkerTrimmed[[x]][[y]])])
					} else {
						LTRed.ends <- end(sampleInfo$sectors[[sector]]$samples[[x]]$LTRed[names(linkerTrimmed[[x]][[y]])])
					}                    
					end(linkerTrimmed[[x]][[y]]) <- end(linkerTrimmed[[x]][[y]]) + LTRed.ends
					start(linkerTrimmed[[x]][[y]]) <- start(linkerTrimmed[[x]][[y]]) + LTRed.ends
					newAttrName <- paste(ifelse(y=="hits","",y),"linkered",sep="")
					sampleInfo$sectors[[sector]]$samples[[x]][[newAttrName]] <- linkerTrimmed[[x]][[y]]
					rm(LTRed.ends)
				}
			}            
			rm("linkerTrimmed")
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
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}}, \code{\link{findPrimers}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' #troubleshootLinkers(sampleInfo,showStats=TRUE)
#'
troubleshootLinkers <- function(sampleInfo, qualityThreshold=0.55, qualityThreshold1=0.75, qualityThreshold2=0.50, doRC=TRUE, parallel=TRUE, samplenames=NULL, ...) {    
    stopifnot(class(sampleInfo)=="SimpleList")
    
    if(!parallel) { registerDoSEQ() }
    
    ## test if there are primed or LTRed sequences in the sampleinfo object ##
    primerFlag <- FALSE
    toProcess <- extractFeature(sampleInfo,feature="LTRed")
    toProcessSamples <- sapply(toProcess,names,simplify=FALSE)
    sectors <- names(toProcess) 
    if(length(sectors)==0) {
        message("No LTRed information found in sampleInfo. Using primed information to find Linkers")
        toProcess <- extractFeature(sampleInfo,feature="primed")
        toProcessSamples <- sapply(toProcess,names,simplify=FALSE)
        sectors <- names(toProcess)
        if(length(sectors)==0) {
            stop("No primed information found in sampleInfo...did you run findPrimers()?")
        }
        primerFlag <- TRUE
    }
    rm(toProcess)
    gc()
    results <- data.frame()
    for(sector in sectors) {
        message("Processing sector ",sector)
        
        ## prepare sample to linker association ##
        sampleLinkers <- toupper(extractFeature(sampleInfo,sector=sector,feature="linkersequence")[[sector]])
        if(length(sampleLinkers)==0 | mean(nchar(sampleLinkers))<=10 | all(is.na(sampleLinkers))) {
            stop("Either Linker sequence is too short (<=10) or no Linkers found in sample information object.")
        }
        
        samplesToProcess <- toProcessSamples[[sector]]
        if(!is.null(samplenames)) {
            samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
        }
        
        toProcess.seqs <- extractSeqs(sampleInfo,sector,samplesToProcess,feature=ifelse(primerFlag,"primed","LTRed"))[[sector]]
        totalSeqs <- sapply(toProcess.seqs,length)
        
        ## do all by all comparison of linkers ##
        message("\tFinding Linkers.")
        for(linkerSeq in unique(as.character(sampleLinkers))) {
            message("Checking ",linkerSeq)
            linkerTrimmed <- foreach(x=iter(samplesToProcess), .inorder=TRUE, .errorhandling="pass", .export=c("linkerSeq","toProcess.seqs","doRC","qualityThreshold","qualityThreshold1","qualityThreshold2","pairwiseAlignSeqs","primerIDAlignSeqs"), .packages="Biostrings") %dopar% {
                if(length(unlist(gregexpr("N",linkerSeq)))>3) {
                    length(primerIDAlignSeqs(toProcess.seqs[[x]], linkerSeq, doRC=doRC, qualityThreshold1=qualityThreshold1, qualityThreshold2=qualityThreshold2, ...)$hits)
                } else {
                    length(pairwiseAlignSeqs(toProcess.seqs[[x]], linkerSeq, "middle", qualityThreshold=qualityThreshold, doRC=doRC, ...))
                }        
            }
            names(linkerTrimmed) <- samplesToProcess
            culprits <- grep("No hits found",linkerTrimmed)
            if(length(culprits)>0) {
                linkerTrimmed[culprits] <- 0
            }   
            results <- rbind(results, data.frame("linkerSeq"=linkerSeq, "samplename"=names(linkerTrimmed), "linkerhits"=as.numeric(unlist(linkerTrimmed)),PercentOfTotal=as.numeric(unlist(linkerTrimmed))/totalSeqs[names(linkerTrimmed)],stringsAsFactors=FALSE))
            cleanit <- gc()
        }        
    }  
    sampleLinkers <- extractFeature(sampleInfo,feature="linkersequence")
    names(sampleLinkers) <- NULL
    sampleLinkers <- unlist(sampleLinkers)
    linkersample <- as.data.frame(sampleLinkers)
    linkersample <- tapply(rownames(linkersample),linkersample$sampleLinkers,paste,collapse=",")

    results$CorrectLinker <- with(results,sampleLinkers[as.character(samplename)]==as.character(linkerSeq))
    results$CorrectSample <- with(results,linkersample[as.character(linkerSeq)])
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
#' @param alignWay method to utilize for detecting the primers. One of following: "slow" (Default), or "fast". Fast, calls \code{\link{vpairwiseAlignSeqs}} and uses \code{\link{vpatternMatch}} at its core, which is less accurate with indels and mismatches but much faster. Slow, calls \code{\link{pairwiseAlignSeqs}} and uses \code{\link{pairwiseAlignment}} at its core, which is accurate with indels and mismatches but slower.
#' @param ... parameters to be passed to \code{\link{pairwiseAlignment}} or \code{\link{vpairwiseAlignSeqs}} depending on which method is defined in 'alignWay' parameter.
#'
#' @return DNAStringSet object with pattern sequence removed from the subject sequences. 
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{extractFeature}}, \code{\link{extractSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findPrimers}}, \code{\link{findLinkers}}
#'
#' @export
#'
#' @examples 
#' #findAndTrimSeq(patternSeq="AGACCCTTTT",subjectSeqs=DNAStringSet(c("AGACCCTTTTGAGCAGCAT","AGACCCTTGGTCGACTCA","AGACCCTTTTGACGAGCTAG")), qualityThreshold=.85, doRC=F, side="left", offBy=1, alignWay = "slow")
#'
findAndTrimSeq <- function(patternSeq, subjectSeqs, side = "left", offBy = 0, alignWay = "slow", ...) {
    
    ## give names to subjectSeqs if not there for matching purpose in trimSeqs()
    removeNamesAfter <- FALSE
    if(is.null(names(subjectSeqs))) {
    	removeNamesAfter <- TRUE
        names(subjectSeqs) <- paste("read",1:length(subjectSeqs))
    }

	coords <- switch(alignWay,
		fast = vpairwiseAlignSeqs(subjectSeqs, patternSeq, side, ...),
		slow = pairwiseAlignSeqs(subjectSeqs, patternSeq, side, ...)
	)
		
	res <- trimSeqs(subjectSeqs, coords, side, offBy)
	if(removeNamesAfter) {
		names(res) <- NULL
	}
	res
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
trimSeqs <- function(dnaSet,coords,side="middle",offBy=0) {
    stopifnot(class(dnaSet) %in% c("DNAStringSet","DNAString"))
    stopifnot(class(coords)=="IRanges")
    
    # check if both dnaSet and coords has 'names' attribute, if yes then check if they have matching names, else check lengths. 
    if(is.null(names(dnaSet)) | is.null(names(coords))) {
        stopifnot(length(dnaSet)==length(coords))
    } else {
        rows <- match(names(coords),names(dnaSet))
        if(any(is.na(rows))) {
            stop("Some of the names in coords are not present in dnaSet")
        }
        if(!is.ordered(rows)) {
            dnaSet <- dnaSet[rows]
            if(!identical(names(dnaSet),names(coords))) {
                stop("Names are not identical between dnaSet and coords parameter")
            }
        }
    }        
    
    # trim by side and check if any of the coords are off the sequence length in dnaSet
    if(tolower(side)=="left") {
        test <- end(coords)+offBy > width(dnaSet) | end(coords)+offBy < 1
        if(any(test)) {
            message("Following sequences were removed from trimming since their coordinates+offBy were out of sequence length: ",paste(names(dnaSet)[test],collapse=", "))
        }
        subseq(dnaSet[!test],start=end(coords[!test])+offBy)
    } else if (tolower(side)=="right") {
        test <- start(coords)-offBy > width(dnaSet) | end(coords)+offBy < 1
        if(any(test)) {
            message("Following sequences were removed from trimming since their coordinates+offBy were out of sequence length: ",paste(names(dnaSet)[test],collapse=", "))
        }
        subseq(dnaSet[!test],end=start(coords[!test])-offBy)
    } else {
        test <- start(coords)+offBy > width(dnaSet) | end(coords)-offBy > width(dnaSet) |  start(coords)+offBy < 1
        if(any(test)) {
            message("Following sequences were removed from trimming since their coordinates+offBy were out of sequence length: ",paste(names(dnaSet)[test],collapse=", "))
        }
        subseq(dnaSet[!test],start=start(coords[!test])+offBy,end=end(coords[!test])-offBy)
    }    
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
getSectorsForSamples <- function(sampleInfo,sector=NULL,samplename=NULL,returnDf=FALSE) {
    stopifnot(class(sampleInfo)=="SimpleList")        
        
    # get all sectors and samplenames in each sector or a specific sector if defined
    if(is.null(sector)) {
        sectors <- names(sampleInfo$sectors)        
    } else {
        sectors <- sector
    }
    samplenames <- sapply(sectors,function(x) names(sampleInfo$sectors[[x]]$samples),simplify=FALSE)
    
    # if specific samplename(s) is defined, then search to find where they are
    if(is.null(samplename)) {
        samplename <- unlist(samplenames)     
    }
    
    if(!all(samplename %in% unlist(samplenames))) {            
        stop("Following sample(s) do not exist on given sector(s) (",paste(sectors,collapse=", "),") in the supplied sampleInfo object: ",paste(samplename[!samplename %in% unlist(samplenames)],collapse=", "))
    }
    sectors <- names(which(unlist(lapply(lapply(samplenames,"%in%",samplename),any)))) 
    if(returnDf) {
        return(do.call(rbind,lapply(sectors, function(x) { 
            data.frame(samplename=samplenames[[x]][samplenames[[x]] %in% samplename], sector=x, stringsAsFactors=FALSE)
        })))
    } else {
        samplenames <- sapply(sectors, function(x) samplenames[[x]][samplenames[[x]] %in% samplename], simplify=FALSE)
        return(list("sectors"=sectors,"samplenames"=samplenames))
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
#'
#' @return a listed DNAStringSet object structed by sector then sample. 
#'
#' @seealso \code{\link{findPrimers}}, \code{\link{findLTRs}}, \code{\link{findLinkers}}, \code{\link{trimSeqs}}, \code{\link{extractFeature}}, \code{\link{getSectorsForSamples}}
#'
#' @export
#'
#' @examples 
#' #extractSeqs(sampleInfo)
#' #extractSeqs(sampleInfo,feature="primed")
#'
extractSeqs <- function(sampleInfo,sector=NULL,samplename=NULL,feature="genomic",trim=TRUE,minReadLength=1,sideReturn=NULL) {
    stopifnot(class(sampleInfo)=="SimpleList")        

    if(is.null(feature)) {
        stop("Please define a feature to extract.")
    }

    # get all sectors and samplenames in each sector or a specific sector
    res <- getSectorsForSamples(sampleInfo,sector,samplename)
    sectors <- res[["sectors"]]
    samplenames <- res[["samplenames"]]
    
    if(feature=="unDecoded") {
       sapply(sectors, function(y) { 
            allmetadata <- metadata(sampleInfo$sectors[[y]])
            if("unDecodedSeqs" %in% names(allmetadata)) {
                allmetadata$unDecodedSeqs 
            } else {
                message("No unDecoded attribute found for the supplied sampleInfo object and sector ",y)
            }
        })
    } else {        
        res <- sapply(sectors, function(y) {
           sapply(samplenames[[y]], function(x,y) { 
                decoded <- sampleInfo$sectors[[y]]$samples[[x]]$decoded
                if(feature!="decoded") {
                    # get ride of ! from feature, else R wont know what to do when making a new object with ! in the front.
                    assign(gsub("!","",feature), sampleInfo$sectors[[y]]$samples[[x]][[gsub("!","",feature)]])
                }
    
                if (feature=="decoded") {
                    decoded
                } else if (feature %in% c("genomic","genomicLinkered")) {
                    primed <- sampleInfo$sectors[[y]]$samples[[x]]$primed
                    LTRed <- sampleInfo$sectors[[y]]$samples[[x]]$LTRed
                    linkered <- sampleInfo$sectors[[y]]$samples[[x]]$linkered
                    
                    if(is.null(LTRed)) {
                        warning("LTRed information not found for",x,"...using primer end as starting boundary.",immediate.=TRUE)                                
                    }

                    if(any(is.null(primed),is.null(linkered))) { 
                        message("No sequences found for requested feature (",feature,") for sample: ",x,"...skipping.") 
                    } else {
                        if(trim) {
                            # get all LTRed reads and make ends = size of each read
                            # if reads dont have LTRs...use primers instead...but throw a warning!
                            if(is.null(LTRed)) {
                                LTRed <- primed
                            }
                            starts <- end(LTRed)
                            starts.name <- names(LTRed)
                            
                            ends <- width(decoded)
                            ends.name <- names(decoded)
                            
                            stopifnot(identical(starts.name,ends.name[ends.name %in% starts.name]))
                            coords <- IRanges(start=starts+1,end=ends[ends.name %in% starts.name],names=starts.name)
                            
                            # alter ends for reads where linker was present.
                            ends <- start(linkered)
                            ends.name <- names(linkered)
                            end(coords[ends.name]) <- ends-1
                            
                            if(feature=="genomicLinkered") { coords <- coords[names(coords) %in% names(linkered)] }
                            
                            # trim it and return non zero length sequences
                            if(length(coords)>0) {
                            	seqs <- trimSeqs(decoded,coords)
								seqs[width(seqs)>=minReadLength]
                            } else {
                            	message("No linkered reads found for sample: ",x,"...skipping.")
                            }                            
                        } else {
                            # just retuning ranges...simply match names from decoded to the request
                                                        
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
                } else {
                    notFeature <- grepl("!",feature)
                    if(notFeature) { ## no need to trim if looking for "not" based feature since there are no coordinates for it
                        toreturn <- decoded[!names(decoded) %in% names(get(gsub("!","",feature)))]
                        if(length(toreturn)) {
                            toreturn
                        } else {
                            message("No sequences found for requested feature (",feature,") for sample: ",x,"...skipping.") 
                        }
                    } else {                    
                        if(is.null(get(feature))) { 
                            message("No sequences found for requested feature (",feature,") for sample: ",x,"...skipping.") 
                        } else {
                            res.seq <- decoded[names(decoded) %in% names(get(feature))]
                            if(trim) {
                                if(is.null(sideReturn)) {
                                    sidetype <- ifelse(grepl("primerID",feature,ignore.case=TRUE), "middle", ifelse(grepl("primed|LTRed",feature,ignore.case=TRUE), "left", ifelse(grepl("linkered",feature,ignore.case=TRUE), "right", "middle")))
                                } else {
                                    sidetype <- tolower(sideReturn)
                                }                            
                                offByLength <- ifelse(sidetype=="middle",0,1)
                                seqs <- trimSeqs(res.seq,get(feature),side=sidetype,offBy=offByLength)
                                seqs[width(seqs)>=minReadLength]
                            } else {
                                res.seq
                            }
                        }
                    }
                }
            }, y=y)
        }, simplify=FALSE)  

        lengthTest <- lapply(lapply(lapply(res, function(x) sapply(x,length)),">",0),which)
        mapply(function(x,y) x[y], res, lengthTest, SIMPLIFY=FALSE)
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
extractFeature <- function(sampleInfo,sector=NULL,samplename=NULL,feature=NULL) {
    stopifnot(class(sampleInfo)=="SimpleList")        

    if(is.null(feature)) {
        stop("Please define a feature to extract.")
    }
    
    # get all sectors and samplenames in each sector or a specific sector
    res <- getSectorsForSamples(sampleInfo,sector,samplename)
    sectors <- res[["sectors"]]
    samplenames <- res[["samplenames"]]
    
    res <- sapply(sectors, function(y) {
           sapply(samplenames[[y]], function(x,y) { 
                if(feature=="metadata") {
                    paste(names(sampleInfo$sectors[[y]]$samples[[x]]), collapse=", ")
                } else {
                    res <- sampleInfo$sectors[[y]]$samples[[x]][[feature]]
                    if(class(res)=="factor") { ## convert any factor based vector to the appropriate regular vector
                        if(!any(is.na(suppressWarnings(as.numeric(levels(res)))))) { as.numeric(as.character(res)) } else { as.character(res) }    
                    } else if (class(res)=="character"){ ## if a numeric vector is stored as character, convert it back to numeric
                        if(!any(is.na(suppressWarnings(as.numeric(res))))) { as.numeric(res) } else { res }    
                    } else {
                        res
                    }
                }
            }, y=y)
    }, simplify=FALSE)
    
    lengthTest <- lapply(lapply(lapply(res, function(x) sapply(x,length)),">",0),which)
    mapply(function(x,y) x[y], res, lengthTest, SIMPLIFY=FALSE)
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
addFeature <- function(sampleInfo,sector=NULL,samplename=NULL,feature=NULL,value=NULL) {
    stopifnot(class(sampleInfo)=="SimpleList")        

    if(is.null(feature)) {
        stop("Please define the new feature(s) to add.")
    }

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
    res <- getSectorsForSamples(sampleInfo,sector,samplename)
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

#' Read fasta/fastq given the path or sampleInfo object.
#'
#' Given a sequence reads file path, the function returns a DNAStringSet object.
#'
#' @param seqFilePath a path to fasta/fastq reads file or a sampleInfo object returned by \code{\link{read.SeqFolder}}
#' @param sector specific sector to reads sequences from. Default is 1, and not required if seqFilePath is a direct file path rather than sampleInfo object.
#'
#' @return a DNAStringSet object.
#'
#' @seealso \code{\link{decodeByBarcode}}, \code{\link{read.SeqFolder}}, \code{\link{extractSeqs}}
#'
#' @export
#'
read.seqsFromSector <- function(seqFilePath=NULL,sector=1) {
    if(is.null(seqFilePath)) {
        stop("Missing seqFilePath!")
    }
    
    if(class(seqFilePath)=="SimpleList") {
        filePath <- normalizePath(grep(paste(sector,seqFilePath$seqfilePattern,sep=""),seqFilePath$seqFilePaths,value=TRUE),mustWork=TRUE)
        if(length(filePath)==0) {
            stop("No sequence file found for sector: ",sector," in seqFilePath variable (",paste(seqFilePath$seqFilePaths,collapse=" * "),") using pattern (",paste(sector,seqFilePath$seqfilePattern,sep=""),")")
        }
        if(length(filePath)>1) {
            stop("Multiple sequence file found for sector: ",sector," in seqFilePath variable (",paste(seqFilePath$seqFilePaths,collapse=" * "),") using pattern (",paste(sector,seqFilePath$seqfilePattern,sep=""),")")
        }
        seqFilePath <- filePath
    }
    
    message("Reading ",seqFilePath)
    if(grepl("fastq$",seqFilePath)) {
    	require(ShortRead)
    	bore <- readFastq(normalizePath(seqFilePath,mustWork=TRUE))
    	dnaSet <- sread(bore)
    	names(dnaSet) <- id(bore)
    	rm(bore)
    } else {
    	dnaSet <- readDNAStringSet(normalizePath(seqFilePath,mustWork=TRUE))
    }
    
    if(length(dnaSet)==0) {
        stop("No sequences found")
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
#' @param format Either fasta (the default) or fastq.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @seealso \code{\link{decodeByBarcode}}, \code{\link{read.SeqFolder}}, \code{\link{extractSeqs}}
#'
#' @note Writing of the files is done using \code{\link{writeXStringSet}} with parameter append=TRUE. If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @export
#'
#' @examples 
#' #write.listedDNAStringSet(dnaSet)
#'
write.listedDNAStringSet <- function(dnaSet, filePath=".", filePrefix="processed", prependSamplenames=TRUE, format="fasta", parallel=TRUE) {
    stopifnot(class(dnaSet)=="list")
    
    if(!parallel) { registerDoSEQ() }
    
    ## do a safety check to see if dnaSet is a list of list or list of DNAStringSet objects ##    
    ## if list of list, then flatten it ##
    if(all(sapply(dnaSet,class)=="list")) {
    	names(dnaSet) <- NULL
    	dnaSet <- unlist(dnaSet)
    }
    
	out <- foreach(i=iter(names(dnaSet)), .inorder=FALSE, .errorhandling="stop", .packages="Biostrings", .export=c("dnaSet","filePath","filePrefix","prependSamplenames")) %dopar% {
		outputSeqs <- dnaSet[[i]]
		if(length(outputSeqs)>0) {
			if(is.null(names(outputSeqs))) {
				message("No names attribute found for ",i," ... using artifically generated names")
				names(outputSeqs) <- paste("read",1:length(outputSeqs),sep="-")
			}
			## remove '.' at the beginning of the filename incase filePrefix is empty
			filename <- gsub("^\\.","",paste(filePrefix,i,ifelse(format=="fastq","fastq","fa"),sep="."))
			filename <- paste(filePath, filename, sep="/")
			if(prependSamplenames) {
				names(outputSeqs) <- paste(i,names(outputSeqs),sep="-")
			}
			writeXStringSet(outputSeqs, file=filename, format=format, append=TRUE) 
		}
	}    
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
    stopifnot(class(sampleInfo)=="SimpleList")
    message("Total sectors:",paste(names(sampleInfo$sectors),collapse=","),"\n")
    do.call(rbind, lapply(names(sampleInfo$sectors), function(sector) {
        res.df <- data.frame(Sector=sector, SampleName=as.character(extractFeature(sampleInfo,sector=sector,feature="samplename")[[sector]]))
        res.df$SampleName <- as.character(res.df$SampleName)
        for (metaD in c("decoded","primed","LTRed","linkered","psl","sites")) {
            res <- sapply(extractFeature(sampleInfo,sector=sector,feature=metaD)[[sector]],length)                        
            if(length(res)>0) {
                res.df[,metaD] <- res[res.df$SampleName]
            } else {
                res.df[,metaD] <- NA
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
summary.elegant <- function(sampleInfo,samplenames=NULL) {
    stopifnot(class(sampleInfo)=="SimpleList")
    cat("Names of all sectors processed:",paste(names(sampleInfo$sectors),collapse=","),"\n")
    for(sector in names(sampleInfo$sectors)) {
        cat(rep("*",50),"\n")
        cat("Sector:",sector,"\n")
        sampleList <- as.character(extractFeature(sampleInfo,sector=sector,feature="samplename")[[sector]])
        if(!is.null(samplenames)) {
            sampleList <- sampleList[sampleList %in% samplenames]
        }
        for(sample.i in sampleList) {
            cat("\t",sample.i,":\n\t")
            metadatalist <- as.character(unlist(strsplit(extractFeature(sampleInfo,sector=sector,samplename=sample.i,feature="metadata")[[sector]],", ")))
            totalDecoded <- extractFeature(sampleInfo,sector=sector,samplename=sample.i,feature="decoded")[[sector]]
            if(length(totalDecoded)>0) {
                totalDecoded <- length(totalDecoded[[1]])
                meta.processed <- grep("ed$",metadatalist,value=TRUE)
                meta.sites <- c("psl","sites")
                for(metadata.i in metadatalist) {
                    cat("| ",metadata.i,": ")
                    feature.i <- extractFeature(sampleInfo,sector=sector,samplename=sample.i,feature=metadata.i)[[sector]][[1]]
                    if(metadata.i %in% meta.processed) {
                        cat(length(feature.i),"(",round(100*(length(feature.i)/totalDecoded),2),"%) ",sep="","\n")
                    } else if(metadata.i %in% meta.sites) {
                        cat(length(unique(feature.i$qName)),"(",round(100*(length(unique(feature.i$qName))/totalDecoded),2),"%) ",sep="","\n")
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
startgfServer <- function(seqDir=NULL, host="localhost", port=5560, gfServerOpts=c(repMatch=112312, stepSize=5, tileSize=10)) {
    if(is.null(seqDir)) {
        stop("Please define the path of nib/2bit files containing the indexed reference sequence(s)")
    }
    
    cmd <- sprintf("gfServer start %s %i %s %s &", 
					host, port, 
					ifelse(!is.null(gfServerOpts),paste(paste("-",names(gfServerOpts),sep=""), gfServerOpts, collapse=" ", sep="="),""), 
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
    
    cmd <- sprintf("kill `ps ax | grep '%s' | grep -v 'grep' | awk '{print $1}'`", paste("gfServer start", host, port))
    system(cmd)
}

#' Align a listed DNAStringSet to a reference using gfClient or standalone BLAT.
#'
#' Align sequences from a listed DNAStringSet object returned from \code{\link{extractSeqs}} to an indexed reference genome using gfServer/gfClient protocol or using standalone BLAT and get a psl file as RangedData object. This function heavily relies on defaults of \code{\link{blatSeqs}}.
#'
#' @param dnaSetList DNAStringSet object containing sequences to be aligned against the reference.
#' @param ... parameters to be passed to \code{\link{blatSeqs}}.
#'
#' @return a list of RangedData object reflecting psl file type per set of sequences.
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{startgfServer}}, \code{\link{stopgfServer}}, \code{\link{blatSeqs}}, \code{\link{read.psl}}, \code{\link{pslToRangedObject}}, \code{\link{read.blast8}}
#'
#' @export
#' 
blatListedSet <- function(dnaSetList=NULL, ...) {
    if(is.null(dnaSetList)) {
        stop("dnaSetList is empty. Please supply a listed DNAStringSet returned from extractSeqs() to be aligned against a reference")
    }
        
    sapply(names(dnaSetList), function(x) {
       do.call(RangedDataList,sapply(names(dnaSetList[[x]]), function(y) {
            if(length(dnaSetList[[x]][[y]])>0) {
                message("BLATing ",y)
                outFiles <- blatSeqs(query=dnaSetList[[x]][[y]], ...)
                read.psl(outFiles, bestScoring=TRUE, asRangedData=TRUE, removeFile=TRUE, parallel=FALSE)
            }
        }))
    })
 }

#' Convert psl dataframe to RangedData/GRanges
#'
#' Convert psl dataframe to RangedData or GRanges object using either the query or target as the reference data column. 
#'
#' @param x dataframe reflecting psl format
#' @param useTargetAsRef use target or query as space or the reference data. Default is TRUE.
#' @param asGRanges make a GRanges object instead of RangedData. Default is FALSE.
#' @param isblast8 the input dataframe blast8 format output from BLAT. Default is FALSE.
#'
#' @return a RangedData/GRanges object reflecting psl file type.
#'
#' @seealso \code{\link{blatListedSet}}
#'
#' @export
#'
#' @examples 
#' #pslToRangedObject(psl)
#' #pslToRangedObject(psl,asGRanges=TRUE)
#' #pslToRangedObject(psl,useTargetAsRef=FALSE)
#'
pslToRangedObject <- function(x, useTargetAsRef=TRUE, asGRanges=FALSE, isblast8=FALSE) {
    if(useTargetAsRef) {
        metadataCols <- c(grep("tName|tStart|tEnd|strand",names(x),invert=TRUE,value=TRUE,fixed=FALSE),ifelse(isblast8,NA,"tStarts"))
        out <- RangedData(space=x$tName,IRanges(start=x$tStart,end=x$tEnd),strand=x$strand, x[,na.omit(metadataCols)])
    } else {
        metadataCols <- c(grep("qName|qStart|qEnd|strand",names(x),invert=TRUE,value=TRUE,fixed=FALSE),ifelse(isblast8,NA,"qStarts"))
        out <- RangedData(space=x$qName,IRanges(start=x$qStart,end=x$qEnd),strand=x$strand, x[,na.omit(metadataCols)])
    }
    
    if(asGRanges) {
    	out <- as(out,"GRanges")
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
splitSeqsToFiles <- function(x, totalFiles=4, suffix="tempy", filename="queryFile.fa") {
    if(is.atomic(x)) {
        message("Splitting file ",x)
        totalSeqs <- length(fasta.info(x,use.names=FALSE))
        chunks <- round(totalSeqs/totalFiles)
        chunks <- ifelse(chunks>0,chunks,totalSeqs) ## incase totalSeqs is lower than number of files to be created!
        starts <- seq(0,totalSeqs,by=chunks) ## create chunks of starts
        for(skippy in starts[starts!=totalSeqs]) {
            filename.out <- paste(x,skippy,suffix,sep=".")
            query.tmp <- readBStringSet(x,nrec=chunks, skip=skippy) ## no need to read the entire file...save memory by reading in N lines
            writeXStringSet(query.tmp,file=filename.out,format="fasta")            
        }
        return(list.files(path=dirname(x), pattern=paste(basename(x),".*",suffix,"$",sep=""), full.names=TRUE))
    } else if (class(x)=="DNAStringSet") {
        message("Splitting Reads.")
        totalSeqs <- length(x)
        chunks <- round(totalSeqs/totalFiles)
        starts <- seq(1,totalSeqs,by=chunks)
        stops <- unique(c(seq(chunks,totalSeqs,by=chunks),totalSeqs))
        stopifnot(length(starts)==length(stops))        
        for(skippy in 1:length(starts)) {
            filename.out <- paste(filename,skippy,suffix,sep=".")            
            writeXStringSet(x[starts[skippy]:stops[skippy]],file=filename.out,format="fasta")            
        }            
        return(list.files(path=".", pattern=paste(filename,".*",suffix,"$",sep=""), full.names=TRUE))
    } else {
        stop("Dont know what is supplied in parameter x.")
    }
}

#' Align sequences using BLAT.
#'
#' Align batch of sequences using standalone BLAT or gfServer/gfClient protocol for alignment against an indexed reference genome. Depending on parameters provided, the function either aligns batch of files to a reference genome using gfClient or takes sequences from query & subject parameters and aligns them using standalone BLAT. If standaloneBlat=FALSE and gfServer is not launched apriori, this function will start one using \code{\link{startgfServer}} and kill it using \code{\link{stopgfServer}} upon successful execution. 
#'
#' @param query an object of DNAStringSet, a character vector, or a path/pattern of fasta files to BLAT. Default is NULL.
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
#' #blatSeqs(dnaSeqs,subjectSeqs,blatParameters=c(minIdentity=90, minScore=10, tileSize=10, dots=10, q="dna", t="dna", out="blast8"))
#' #blatSeqs(dnaSeqs,"/usr/local/blatSuite34/hg18.2bit",standaloneBlat=FALSE)
#' #blatSeqs("mySeqs.fa","/usr/local/blatSuite34/hg18.2bit",standaloneBlat=FALSE)
#'
blatSeqs <- function(query=NULL, subject=NULL, standaloneBlat=TRUE, port=5560, host="localhost", parallel=TRUE, gzipResults=TRUE, blatParameters=c(minIdentity=70, minScore=5, stepSize=5, tileSize=10, repMatch=112312, dots=50, q="dna", t="dna", out="psl")) {

    ## get all BLAT options from the system for comparison to blatParameters later
    suppressWarnings(blatOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", grep("\\s+-.+=", system("blat",intern=TRUE), value=TRUE))))
    suppressWarnings(gfClientOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", grep("\\s+-.+=", system("gfClient",intern=TRUE), value=TRUE))))

    if(!standaloneBlat) { 
        message("Using gfClient protocol to perform BLAT.")
        if(is.null(port)) {
            stop("The port paramter is empty. Please define the port used to start gfServer with")
        }
    }

    ## check the subject parameter
    if(is.null(subject)) {
        stop("The subject parameter is empty. Please supply subject sequences or a path to 2bit or nib files to serve as reference/target")
    } else {
        subjectFile <- NULL
        if(is.atomic(subject)) {
            if (any(grepl("\\.2bit$|\\.nib$",subject))) {
                if(standaloneBlat) { stop("Standalone BLAT cannot be used when subject is an indexed nib or 2bit file.") }
                indexFileDir <- dirname(subject)
                subjectFile <- list.files(path=indexFileDir, pattern=basename(subject), full.names=TRUE)
                if(length(subjectFile)==0) { stop("The file(s) supplied in subject parameter doesn't exist.") }
            } else {
                ## change object type if necessary for troubleshooting purpose in later steps
                subject <- DNAStringSet(subject)
            }
        }
        
        if(is.null(subjectFile)) {
            ## subjectFile is still null so it means that subject is a DNAStringSet
            if(is.null(names(subject))) { ## add names of subject if not present
                names(subject) <- paste("subject",1:length(subject))
            }

            ## write out the subject sequences into a fasta file
            filename.seq <- "subjectFile.fa.tempyS"
            writeXStringSet(subject,file=filename.seq,format="fasta")                                  
            subjectFile <- filename.seq
        }
    }
            
    ## check the query parameter
    if(is.null(query)) {
        stop("The query parameter is empty. Please supply reads to be aligned")
    } else {
        queryFiles <- NULL
        if(is.atomic(query)) {
           if (any(grepl("\\.fna$|\\.fa$|\\*",query))) {
                queryFiles <- list.files(path=dirname(query), pattern=basename(query), full.names=TRUE)            
                if(parallel) {
                    ## split the fasta files into smaller chunks for parallel BLATing
                    queryFiles <- unlist(sapply(queryFiles,function(f) splitSeqsToFiles(f,getDoParWorkers(),"tempyQ")), use.names=FALSE)                    
                }
            } else {
                ## change object type if necessary for troubleshooting purpose in later steps
                query <- DNAStringSet(query)
            }
        }
        
        if(is.null(queryFiles)) {
            ## queryFiles is still null so it means that query is a DNAStringSet           
            if(is.null(names(query))) {  ## fix names of query if not present
                names(query) <- paste("read",1:length(query))
            }  
            
            ## write out the query sequences into fasta files
            if(parallel) {
                queryFiles <- splitSeqsToFiles(query,getDoParWorkers(),"tempyQ")
            } else {
                queryFiles <- "queryFile.fa.tempyQ"
                writeXStringSet(query,file=queryFiles,format="fasta")                
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
        filenames <- foreach(x=iter(queryFiles),.inorder=FALSE,.export=c("blatOpts","subjectFile","gzipResults")) %dopar% {
            filename.out <- paste(x,blatOpts["out"],sep=".")
            cmd <- paste("blat", paste(paste("-",names(blatOpts),sep=""), blatOpts, collapse=" ", sep="="), "-noHead", subjectFile, x, filename.out)
            message(cmd)
            system(cmd)
            if(grepl("\\.tempyQ$",x)) { system(paste("rm",x)) } ## no need to save splitted files!
            if(gzipResults) { system(paste("gzip",filename.out)); filename.out <- paste(filename.out,"gz",sep=".") }
            filename.out
        }
        if(grepl("\\.tempyS$",subjectFile)) { system(paste("rm",subjectFile)) }
    } else {
		# start the gfServer if not started already! #
		killFlag <- FALSE
		searchCMD <- sprintf("gfServer status %s %s", host, port)
        if(system(searchCMD,ignore.stderr=TRUE)!=0) {
        	message("Starting gfServer.")        
            startgfServer(seqDir=subjectFile, host=host, port=port, 
            			  gfServerOpts=c(repMatch=blatParameters[['repMatch']], 
            			                 stepSize=blatParameters[['stepSize']], 
            			                 tileSize=blatParameters[['tileSize']]))
            killFlag <- TRUE
        }         

        gfClientOpts <- blatParameters[names(blatParameters) %in% gfClientOpts]
        stopifnot(length(subjectFile)>0)
        filenames <- foreach(x=iter(queryFiles), .inorder=FALSE, .export=c("gfClientOpts","host","port","indexFileDir","gzipResults")) %dopar% {
            filename.out <- paste(x,gfClientOpts["out"],sep=".")
            cmd <- paste("gfClient",
            			 paste(paste("-",names(gfClientOpts),sep=""), gfClientOpts, collapse=" ", sep="="), 
            			 "-nohead", host, port, "/", x, filename.out)
            message(cmd)
            system(cmd)
            if(grepl("\\.tempyQ$",x)) { system(paste("rm",x)) } ## no need to save splitted files!
            if(gzipResults) { system(paste("gzip",filename.out)); filename.out <- paste(filename.out,"gz",sep=".") }
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
#' Given filename(s), the function reads the psl file format from BLAT as a data frame and performs basic score filtering if indicated. Any other file format will yield errors or erroneous results.
#'
#' @param pslFile psl filename, or vector of filenames, or a pattern of files to import.
#' @param bestScoring report only best scoring hits instead of all hits. Default is TRUE. Score is calculated by matches-misMatches-qBaseInsert-tBaseInsert.
#' @param asRangedData return a RangedData object instead of a dataframe. Default is TRUE. Saves memory!
#' @param removeFile remove the psl file(s) after importing. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @return a dataframe or RangedData object reflecting psl file type.
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{startgfServer}}, \code{\link{blatSeqs}}, \code{\link{read.blast8}}
#'
#' @export
#'
#' @examples 
#' #read.psl(pslFile="processed.*.psl$")
#' #read.psl(pslFile=c("sample1hits.psl","sample2hits.psl"))
#'
read.psl <- function(pslFile=NULL, bestScoring=TRUE, asRangedData=FALSE, removeFile=TRUE, parallel=FALSE) {
    if(is.null(pslFile)) {
        stop("pslFile parameter empty. Please supply a filename to be read.")
    }
    
    if (any(grepl("\\*",pslFile))) {
        pslFile <- list.files(path=dirname(pslFile), pattern=basename(pslFile), full.names=TRUE)  ## vector of filenames          
    }
    
    if(length(pslFile)==0) { stop("No file(s) found with given paramter in pslFile:",pslFile) }
    
    if(!parallel) { registerDoSEQ() }
    
    ## setup psl columns + classes
    cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts")
    cols.class <- c(rep("numeric",8),rep("character",2),rep("numeric",3),"character",rep("numeric",4),rep("character",3))
    
    hits <- foreach(x=iter(pslFile), .inorder=FALSE, .export=c("cols","cols.class","bestScoring","removeFile")) %dopar% {
        message(x)
        hits.temp <- read.delim(x, header=FALSE, col.names=cols, stringsAsFactors=FALSE, colClasses=cols.class)
        if(removeFile) { system(paste("rm", x)) }
        if(bestScoring) {  ## do round one of bestScore here to reduce file size          
            hits.temp$score <- with(hits.temp,matches-misMatches-qBaseInsert-tBaseInsert)
            hits.temp <- hits.temp[with(hits.temp, order(qName, -score)), ]
            bestScore <- with(hits.temp,tapply(score,qName,max))
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

    message("Ordering by qName and cherry picking!")
    hits <- arrange(hits,qName)
    
    if(bestScoring) { ## do round two of bestScore incase any got missed in round one                
        hits$score <- with(hits,matches-misMatches-qBaseInsert-tBaseInsert)
        hits <- hits[with(hits, order(qName, -score)), ]
        bestScore <- with(hits,tapply(score,qName,max))
        isBest <- with(hits, score==bestScore[qName])
        hits <- hits[isBest,]
        rm("isBest","bestScore")
    }
    
    if(asRangedData) {
        hits <- pslToRangedObject(hits, useTargetAsRef=TRUE)
    }
    
    return(hits)
}

#' Read blast8 file(s) outputted by BLAT
#'
#' Given filename(s), the function reads the blast8 file format from BLAT as a data frame and performs basic score filtering if indicated. Any other file format will yield errors or erroneous results.
#'
#' @param files blast8 filename, or vector of filenames, or a pattern of files to import.
#' @param asRangedData return a RangedData object instead of a dataframe. Default is TRUE. Saves memory!
#' @param removeFile remove the blast8 file(s) after importing. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @return a dataframe or RangedData object reflecting blast8 file type.
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{startgfServer}}, \code{\link{blatSeqs}}
#'
#' @export
#'
#' @examples 
#' #read.blast8(files="processed.*.blast8$")
#' #read.blast8(files=c("sample1hits.blast8","sample2hits.blast8"))
#'
read.blast8 <- function(files=NULL, asRangedData=FALSE, removeFile=TRUE, parallel=FALSE) {
    if(is.null(files)) {
        stop("files parameter empty. Please supply a filename to be read.")
    }
    
    if (any(grepl("\\*",files))) {
        files <- list.files(path=dirname(files), pattern=basename(files), full.names=TRUE)  ## vector of filenames          
    }
    
    if(length(files)==0) { stop("No file(s) found with given paramter in files:",files) }
    
    if(!parallel) { registerDoSEQ() }
    
    ## setup blast8 columns + classes
    cols <- c("qName", "tName", "identity", "span", "misMatches", "gaps", "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")
	cols.class<-c(rep("character",2),rep("numeric",10))
    
    hits <- foreach(x=iter(files), .inorder=FALSE, .export=c("cols","cols.class","bestScoring","removeFile")) %dopar% {
        message(x)
        hits.temp <- read.delim(x, header=FALSE, col.names=cols, stringsAsFactors=FALSE, colClasses=cols.class)
        hits.temp$strand <- with(hits.temp,ifelse(tStart>tEnd,"-","+"))
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

    message("Ordering by qName and cherry picking!")
    hits <- arrange(hits,qName)
    
    if(asRangedData) {
        hits <- pslToRangedObject(hits, useTargetAsRef=TRUE, isblast8=TRUE)
    }
    
    return(hits)
}

#' Obtain integration sites from BLAT output
#'
#' Given a RangedData object from \code{\link{read.psl}}, the function uses specified filtering parameters to obtain integration sites and maintain sequence attrition. The function will remove any non-best scoring alignments from the object if not already filtered apriori.
#'
#' @param psl.rd a RangedData object reflecting psl format where tName is the spaces.
#' @param startWithin upper bound limit of where the alignment should start within the query. Default is 3.
#' @param alignRatioThreshold cuttoff for (alignment span/read length). Default is 0.7.
#' @param genomicpercentidentity cuttoff for (1-(misMatches/matches)). Default is 0.98.
#' @param correctByqStart use qStart to correct genomic position. This would account for sequencing/trimming errors. Position=ifelse(strand=="+",tStart-qStart,tEnd+qStart). Default is TRUE.
#'
#' @return a RangedData object with integration sites which passed all filtering criteria. Each filtering parameter creates a new column to flag if a sequence/read passed that filter which follows the scheme: 'pass.FilterName'.
#'
#' @seealso \code{\link{startgfServer}}, \code{\link{read.psl}}, \code{\link{blatSeqs}}, \code{\link{blatListedSet}}, \code{\link{findIntegrations}}, \code{\link{pslToRangedObject}}, \code{\link{clusterSites}}, \code{\link{otuSites2}}, \code{\link{crossOverCheck}}, \code{\link{read.blast8}}
#'
#' @export
#'
#' @examples 
#' #getIntegrationSites(test.psl.rd)
#'
getIntegrationSites <- function(psl.rd=NULL, startWithin=3, alignRatioThreshold=0.7, genomicpercentidentity=0.98, correctByqStart=TRUE) {
    stopifnot(class(psl.rd)=="RangedData" & !is.null(psl.rd) & !is.null(startWithin) & !is.null(alignRatioThreshold) & !is.null(genomicpercentidentity))

    ## get the integration position by correcting for any insertions due to sequencing errors ##    
    if(correctByqStart) {
        psl.rd$Position <- ifelse(psl.rd$strand=="+",start(psl.rd)-psl.rd$qStart,end(psl.rd)+psl.rd$qStart)
    } else {
        psl.rd$Position <- ifelse(psl.rd$strand=="+",start(psl.rd),end(psl.rd))
    }
    
    # get scores for picking best hits and identify multihits later
    # check if scoring filtering hasn't already been applied by blat functions
    if(!"score" %in% colnames(psl.rd)) {
        psl.rd$score <- with(psl.rd,matches-misMatches-qBaseInsert-tBaseInsert)
        bestScore <- tapply(psl.rd$score,space(psl.rd),max)
        isBest <- with(psl.rd, score==bestScore[space(psl.rd)])
        psl.rd <- psl.rd[isBest,]
        rm("isBest","bestScore")
        gc()
    }    
    
    message("Performing QC checks.")
    # remove rows where the best hit dont start within first X bp
    psl.rd$pass.startWithin <- psl.rd$qStart<=startWithin

    # check if aligned ratio matches the threshold
    psl.rd$alignRatio <- with(psl.rd,score/qSize)
    psl.rd$pass.alignRatio <- psl.rd$alignRatio >= alignRatioThreshold
    
    # check for %identity    
    psl.rd$percIdentity <- with(psl.rd,1-(misMatches/matches))  
    psl.rd$pass.percIdentity <- psl.rd$percIdentity >= genomicpercentidentity

    ## find which query aligned to multiple places with equally good score aka...multihits
    cloneHits <- table(psl.rd$qName)
    psl.rd$isMultiHit <- as.logical(cloneHits[as.character(psl.rd$qName)]>1)
    rm(cloneHits)    
    gc()
    
    psl.rd$pass.allQC <- with(psl.rd,pass.percIdentity & pass.alignRatio & pass.startWithin)
        
    return(psl.rd)
}

#' Cluster/Correct values within a window based on their frequency given discrete factors
#'
#' Given a group of discrete factors (i.e. position ids) and integer values, the function tries to correct/cluster the integer values based on their frequency in a defined windowsize.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand). Required if psl.rd parameter is not defined. 
#' @param value a vector of integer with values that needs to corrected/clustered (i.e. Positions). Required if psl.rd parameter is not defined. 
#' @param grouping additional vector of grouping by which to pool the rows (i.e. samplenames). Default is NULL.
#' @param psl.rd a RangedData object returned from \code{\link{getIntegrationSites}}. Default is NULL. 
#' @param weight a numeric vector of weights to use when calculating frequency of value by posID and grouping if specified. Default is NULL.
#' @param windowSize size of window for which values should be corrected/clustered. Default is 5.
#' @param byQuartile flag denoting whether quartile based technique should be employed. See notes for details. Default is TRUE.
#' @param quartile if byQuartile=TRUE, then the quartile which serves as the threshold. Default is 0.70.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}. Process is split by the grouping the column.
#'
#' @note The algorithm for clustering when byQuartile=TRUE is as follows: for all values in each grouping, get a distribution and test if their frequency is >= quartile threshold. For values below the quartile threshold, test if any values overlap with the ones that passed the threshold and is within the defined windowSize. If there is a match, then merge with higher value, else leave it as is. This is only useful if the distribution is wide and polynodal. When byQuartile=FALSE, for each group the values within the defined window are merged with the next highest frequently occuring value, if freuquencies are tied then lowest value is used to represent the cluster. When psl.rd is passed, then multihits are ignored and only unique sites are clustered. All multihits will be tagged as a good 'clusterTopHit'.
#'
#' @return a data frame with clusteredValues and frequency shown alongside with the original input. If psl.rd parameter is defined then a RangedData object is returned with three new columns appended at the end: clusteredPosition, clonecount, and clusterTopHit (a representative for a given cluster chosen by best scoring hit!). 
#'
#' @seealso \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, \code{\link{otuSites}}, \code{\link{otuSites2}}, \code{\link{crossOverCheck}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #clusterSites(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), value=c(rep(1000,2),5832,1000,12324,65738,928042), grouping=c('a','a','a','b','b','b','c'))
#' #clusterSites(grouping=test.psl.rd$grouping, psl.rd=test.psl.rd)
#'
clusterSites <- function(posID=NULL, value=NULL, grouping=NULL, psl.rd=NULL, weight=NULL, windowSize=5, byQuartile=FALSE, quartile=0.70, parallel=TRUE) {
	if(!parallel) { registerDoSEQ() }      
    if(is.null(psl.rd)) {
        stopifnot(!is.null(posID))
        stopifnot(!is.null(value))
    } else {
        ## dereplicate & converge sites! ##
        if(!"Position" %in% colnames(psl.rd)) {
            stop("The object supplied in psl.rd parameter does not have Position column in it. Did you run getIntegrationSites() on it?")
        }
		
		good.row <- rep(TRUE,nrow(psl.rd))
		isthere <- grepl("isMultiHit",colnames(psl.rd),ignore.case=TRUE)		
		if(any(isthere)) { ## see if multihit column exists
			message("Found 'isMultiHit' column in the data. These rows will be ignored for the calculation.")
			good.row <- good.row & !psl.rd[[which(isthere)]]
		}

        posIDs <- paste0(space(psl.rd),psl.rd$strand)
        values <- psl.rd$Position
        if(is.null(weight)) { ## see if sequences were dereplicated before in the pipeline which adds counts=x identifier to the deflines
            weight <- suppressWarnings(as.numeric(sub(".+counts=(\\d+)","\\1",psl.rd$qName)))
            if(all(is.na(weight))) { weight <- NULL } else { weight[is.na(weight)] <- 1 }
        }
        grouping <- if(is.null(grouping)) { "" } else { grouping }
        clusters <- clusterSites(posIDs[good.row], 
        						 values[good.row], 
        						 grouping=grouping[good.row], 
        						 weight=weight, 
								 windowSize=windowSize, 
        						 byQuartile=byQuartile, quartile=quartile)
        
        message("Adding clustered data back to psl.rd.")        
        clusteredValues <- with(clusters,split(clusteredValue,paste0(posID,value,grouping)))
        groupingVals <- paste0(posIDs, values, grouping)[good.row]
        psl.rd$clusteredPosition <- psl.rd$Position
        psl.rd$clusteredPosition[good.row] <- as.numeric(clusteredValues[groupingVals])
        
        ## add frequency of new clusteredPosition ##
        clusteredValueFreq <- with(clusters, split(clusteredValue.freq, paste0(posID,value,grouping)))
        psl.rd$clonecount <- 0
        psl.rd$clonecount[good.row] <- as.numeric(clusteredValueFreq[groupingVals])
        rm("clusteredValueFreq","clusteredValues","clusters")
        cleanit <- gc()
        
        ## pick best scoring hit to represent a cluster ##
        message("Picking best scoring hit to represent a cluster.")
        isthere <- grepl("score",colnames(psl.rd),ignore.case=TRUE)		
        if(!any(isthere)) {
            message("No 'score' column found in the data. Using 'qSize' as an alternative.")
            isthere <- grepl("qSize",colnames(psl.rd),ignore.case=TRUE)		
            if(!any(isthere)) {
            	stop("No 'qSize' column found in the data either...can't pick the best hit :(")
            }            
        }
        groupingVals <- paste0(space(psl.rd), psl.rd$strand, psl.rd$clusteredPosition, grouping)
        bestScore <- tapply(psl.rd[[which(isthere)]], groupingVals, max)
        isBest <- psl.rd[[which(isthere)]]==bestScore[groupingVals]
        
        ## pick the first match for cases where >1 reads with the same coordinate had the same best scores ##
        tocheck <- which(isBest)
        res <- tapply(tocheck,names(tocheck),"[[",1) 
        
        psl.rd$clusterTopHit <- FALSE
        psl.rd$clusterTopHit[res] <- TRUE
        psl.rd$clusterTopHit[!good.row] <- TRUE
        
        message("Cleaning up!")
        rm("isBest","bestScore","posIDs","values","groupingVals")
        cleanit <- gc()
        return(psl.rd)
    }
    
    # get frequencies of each posID & value combination by grouping #
    groups <- if(is.null(grouping)) { "" } else { grouping }
    weight2 <- if(is.null(weight)) { 1 } else { weight }
    sites <- arrange(data.frame(posID, value, grouping=groups, weight=weight2, posID2=paste0(groups,posID), stringsAsFactors=FALSE), posID2, value)
    sites <- count(sites, c("posID","value","grouping","posID2"), wt_var="weight")    
    rm("groups","weight2")
    
    if(byQuartile) {
        message("Clustering by quartile: ",quartile)
        # obtain the defined quartile of frequency per posID & grouping #
        sites <- arrange(sites, posID2, value, desc(freq))
        quartiles <- with(sites,tapply(freq,posID2,quantile,probs=quartile,names=FALSE))
        sites$belowQuartile <- with(sites,freq < quartiles[posID2])
        rm(quartiles)
        
        if(any(sites$belowQuartile)) {
            # for values belowQuartile, see if any within defined windowSize of aboveQuartile #
            pos.be <- with(subset(sites,belowQuartile,drop=TRUE),GRanges(IRanges(start=value,width=1),seqnames=posID2, freq=freq))
            pos.ab <- with(subset(sites,!belowQuartile,drop=TRUE),GRanges(IRanges(start=value,width=1),seqnames=posID2, freq=freq))
            pos.overlap <- as.data.frame(as.matrix(findOverlaps(pos.be,pos.ab,maxgap=windowSize,ignore.strand=TRUE)))
                
            # for overlapping values, merge them with the biggest respective aboveQuartile site #
            pos.overlap$freq <- values(pos.ab[pos.overlap[,"subjectHits"]])$freq
            maxs <- with(pos.overlap,tapply(freq,as.character(queryHits),max))
            pos.overlap$isMax <- with(pos.overlap, freq == maxs[as.character(queryHits)])
            rm(maxs)
            
            pos.overlap <- subset(pos.overlap,isMax, drop=TRUE)
            
            # if there are >1 biggest respective aboveQuartile site, then choose the closest one
            # if tied, then use the latter to represent the site
            counts <- xtabs(isMax~queryHits,pos.overlap)
            if(length(table(counts))>1) {
                toFix <- as.numeric(names(which(counts>1)))
                rows <- pos.overlap$queryHits %in% toFix
                pos.overlap$aboveQuartileValue <- pos.overlap$belowQuartileValue <- pos.overlap$valueDiff <- 0
                pos.overlap$aboveQuartileValue[rows] <- start(pos.ab[pos.overlap[rows,"subjectHits"]])
                pos.overlap$belowQuartileValue[rows] <- start(pos.be[pos.overlap[rows,"queryHits"]])
                pos.overlap$valueDiff[rows] <- with(pos.overlap[rows,],abs(aboveQuartileValue-belowQuartileValue))
                mins <- with(pos.overlap[rows,],tapply(valueDiff,as.character(queryHits),min))
                pos.overlap$isClosest <- TRUE
                pos.overlap$isClosest[rows] <- with(pos.overlap[rows,], valueDiff == mins[as.character(queryHits)])
                pos.overlap <- subset(pos.overlap,isMax & isClosest,drop=TRUE)
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
                res <- clusterSites(rep(as.character(seqnames(loners)),times=times.rep),rep(start(loners),times=times.rep),byQuartile=FALSE)
            }
            pos.overlap <- rbind(pos.overlap[,c("queryHits","clusteredValue")], 
                                 data.frame(queryHits=rows, clusteredValue=as.numeric(res$clusteredValue)))
        
            sites$clusteredValue <- sites$value
            sites$clusteredValue[sites$belowQuartile][pos.overlap[,"queryHits"]] <- pos.overlap$clusteredValue
            stopifnot(any(!is.na(sites$clusteredValue)))            
        } else {
            message("No sites found below defined quartile. Try to increase the quartile or use standard clustering, byQuartile=FALSE.")
        }
    } else {
        message("Clustering by minimum overlap.")
        
        sites <- split(sites, sites$grouping)
        
        sites <- foreach(x=iter(sites), .inorder=FALSE, .packages="IRanges", .combine=rbind) %dopar% {
        	
        	## find overlapping positions using findOverlaps() using maxgap adjusted by windowSize! ##
			sites.rl <- with(x,RangedData(space=posID2,IRanges(start=value,width=1),freq))
			
			# the key part is ignoreSelf=TRUE,ignoreRedundant=FALSE...helps overwrite values at later step
			res <- as.data.frame(as.matrix(findOverlaps(sites.rl,ignoreSelf=TRUE,ignoreRedundant=FALSE,select="all",maxgap=windowSize))) 
			
			# add accessory columns to dictate decision making!
			# q = query, s = subject, val = value, freq = frequency of query/subject
			res$q.val <- start(sites.rl)[res$queryHits]; res$s.val <- start(sites.rl)[res$subjectHits]
			res$q.freq <- sites.rl$freq[res$queryHits]; res$s.freq <- sites.rl$freq[res$subjectHits]
			res$dist <- with(res,abs(q.val-s.val))
			stopifnot(!any(res$dist>5)) ## do safety checking!
	
			# favor a lower value where frequence/cloneCount is tied, else use the value of the highest frequency!
			res$val <- with(res,ifelse(q.freq==s.freq, ifelse(q.val < s.val,q.val,s.val), ifelse(q.freq >= s.freq,q.val,s.val))) 
			
			# for cases where there are >1 matches between query & subject...find the one with the highest frequency and merge with that.
			# if all frequencies are the same, then use the lowest value to represent the cluster!
			res$maxFreq <- with(res,pmax(q.freq,s.freq))    
			maxes <- with(res,tapply(maxFreq,queryHits,max))
			res$ismaxFreq <- with(res,maxFreq==maxes[as.character(queryHits)])        
			res <- arrange(res,desc(queryHits),desc(val)) ## VIP step...this is what merges high value to low value for ties in the hash structure below!!!
			hash.df <- unique(subset(res,ismaxFreq)[,c("queryHits","val")])
			clustered <- structure(as.numeric(hash.df$val),names=as.character(hash.df$queryHits))
			rm(hash.df)
			
			# trickle results back to sites
			x$clusteredValue <- x$value
			x$clusteredValue[as.numeric(names(clustered))] <- as.numeric(clustered)
			rm("clustered","res")
			cleanit <- gc()
			x
        }        
    }
    
    message("\t - Adding clustered value frequencies.")
    # get frequency of clusteredValue
    counts <- count(sites[,-grep("value",names(sites),fixed=TRUE)],c("posID2","clusteredValue"),wt_var="freq")
    names(counts)[grep("freq",names(counts),fixed=TRUE)] <- "clusteredValue.freq"
    sites <- merge(sites,counts)
    
    if(byQuartile) {
        sites <- sites[,c("posID","value","freq","clusteredValue","clusteredValue.freq","grouping")]
    }
    
    sites$posID2<-NULL
    if(is.null(grouping)) { sites$grouping<-NULL }
    if(is.null(weight)) { sites$weight<-NULL }
    
    return(sites)
}

#' Make OTUs of discrete positions grouped by reads
#'
#' Given a group of discrete positions per read/clone, the function tries to yield a unique OTU ID for the collection based on overlap of discrete positions to other reads/clones by grouping. This is mainly useful when each readID has many posID which needs to be considered as one single group of sites.
#'
#' @param posID a vector of discrete positions, i.e. Chr,strand,Position.
#' @param readID a vector of read/clone names which is unique to each row, i.e. deflines.
#' @param grouping additional vector of grouping by which to pool the rows (i.e. samplenames). Default is NULL.
#' @param psl.rd a RangedData object returned from \code{\link{clusterSites}}. Default is NULL. 
#'
#' @note The algorithm for making OTUs of sites is as follows: for each readID check how many positions are there. Separate readIDs with only position from the rest. Check if any readIDs with >1 position match to any readIDs with only one position. If there is a match, then assign both readIDs with the same OTU ID. Check if any positions from readIDs with >1 position match any other readIDs with >1 position. If yes, then assign same OTU ID to all readIDs sharing 1 or more positions. 
#'
#' @return a data frame with posID, readID, grouping, and otuID. If psl.rd parameter is defined, then a RangedData object where object is first filtered by clusterTopHit column and the otuID column appended at the end.
#'
#' @seealso \code{\link{clusterSites}}, \code{\link{otuSites2}}, \code{\link{crossOverCheck}}, \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #otuSites(posID=c('chr1-1000','chr1-1000','chr2-1000','chr2+1000','chr15-1000','chr16-1000','chr11-1000'), readID=paste('read',sample(letters,7),sep='-'), grouping=c('a','a','a','b','b','b','c'))
#' #otuSites(psl.rd=test.psl.rd)
#'
otuSites <- function(posID=NULL, readID=NULL, grouping=NULL, psl.rd=NULL, parallel=TRUE) {
  if(!parallel) { registerDoSEQ() }
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(readID))
  } else {
    ## find the otuID by clusters ##
    if(!"clusterTopHit" %in% colnames(psl.rd) | !"clusteredPosition" %in% colnames(psl.rd)) {
      stop("The object supplied in psl.rd parameter does not have 'clusterTopHit' or 'clusteredPosition' column in it. 
            	  Did you run clusterSites() on it?")
    }
    
    good.rows <- psl.rd$clusterTopHit
    posIDs <- paste(space(psl.rd),psl.rd$strand,psl.rd$clusteredPosition,sep="")
    if("qName" %in% colnames(psl.rd)) {
    	readID <- psl.rd$qName
    } else if ("Sequence" %in% colnames(psl.rd)) {
    	readID <- psl.rd$Sequence
    } else {
    	stop("No readID type column found in psl.rd object.")
    }
    grouping <- if(is.null(grouping)) { rep("A",nrow(psl.rd)) } else { grouping }
    
    otus <- otuSites(posIDs[good.rows], readIDs[good.rows], grouping[good.rows])
    
    message("Adding otuIDs back to psl.rd.")        
    otuIDs <- with(otus,split(otuID,paste(posID,readID,grouping)))
    psl.rd$otuIDs <- NA
    psl.rd$otuIDs[good.rows] <- as.numeric(otuIDs[paste(posIDs, readIDs, 
    													grouping)[good.rows]])
    
    rm("otus","otuIDs","posIDs","readIDs","grouping")
    cleanit <- gc()
    return(psl.rd)
  }
  
  groups <- if(is.null(grouping)) { "A" } else { grouping }
  sites <- data.frame(posID,readID,grouping=groups,stringsAsFactors=FALSE)
  rm(groups)
  
  ## get unique posIDs per readID by grouping 
  reads <- ddply(sites, .(grouping,readID), summarise, posIDs=paste(sort(unique(posID)),collapse=","), counts=length(unique(posID)), .parallel=parallel)
  reads <- arrange(reads, grouping, posIDs)
  
  # create initial otuID by assigning a numeric ID to each collection of posIDs per grouping
  reads$otuID <- unlist(lapply(lapply(with(reads,split(posIDs,grouping)),as.factor),as.numeric)) 
  reads$newotuID <- reads$otuID 
  reads$check <- TRUE
  
  ## see if readID with a unique or single posID matches up to a readID with >1 posIDs, if yes then merge
  singles <- reads$counts==1
  if(any(singles)) {
    message('Merging non-singletons with singletons.')
    toCheck <- with(reads[singles,], split(posIDs,grouping))
    toCheck.ids <- with(reads[singles,], split(otuID,grouping))
    toCheck <- sapply(names(toCheck), function(x) { names(toCheck[[x]]) <- toCheck.ids[[x]]; toCheck[[x]] } )
    rm(toCheck.ids)
    
    allposIDs <- with(reads[!singles,],split(posIDs,grouping))
    
    for(f in intersect(names(toCheck),names(allposIDs))) {
      query <- structure(paste(toCheck[[f]],",",sep=""), names=names(toCheck[[f]])) # this is crucial to avoid matching things like xyzABC to xyz
      subject <- paste(allposIDs[[f]],",",sep="") # this is crucial to avoid matching things like xyzABC to xyz
      res <- sapply(query, grep, x=subject, fixed=TRUE)
      res <- res[sapply(res,length)>0]                    
      res <- structure(unlist(res,use.names=F),names=rep(names(res),sapply(res,length)))
      reads[!singles & reads$grouping==f,"newotuID"][as.numeric(res)] <- as.numeric(names(res))
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
        res <- which(unlist(lapply(lapply(strsplit(tocheck$posIDs,","),"%in%", unlist(strsplit(posId,","))),any)))
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
  ots.ids <- with(reads,split(newotuID,paste0(readID,grouping)))
  if(!is.numeric(ots.ids)) {
    stop("Something went wrong merging non-singletons. Multiple OTUs assigned to one readID most likely!")
  }
  sites$otuID <- as.numeric(unlist(ots.ids[with(sites,paste0(readID,grouping))]))
  
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
#' @param psl.rd a RangedData object returned from \code{\link{clusterSites}}. Default is NULL. 
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}. Process is split by the grouping the column.
#'
#' @note The algorithm for making OTUs of sites is as follows: for each readID check how many positions are there. Separate readIDs with only position from the rest. Check if any readIDs with >1 position match to any readIDs with only one position. If there is a match, then assign both readIDs with the same OTU ID. Check if any positions from readIDs with >1 position match any other readIDs with >1 position. If yes, then assign same OTU ID to all readIDs sharing 1 or more positions.
#'
#' @return a data frame with binned values and otuID shown alongside the original input. If psl.rd parameter is defined, then a RangedData object where object is first filtered by clusterTopHit column and the otuID column appended at the end.
#'
#' @seealso \code{\link{clusterSites}}, \code{\link{otuSites}}, \code{\link{crossOverCheck}}, \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #otuSites2(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), value=c(rep(1000,2),5832,1000,12324,65738,928042), readID=paste('read',sample(letters,7),sep='-'), grouping=c('a','a','a','b','b','b','c'))
#' #otuSites2(psl.rd=test.psl.rd)
#'
otuSites2 <- function(posID=NULL, value=NULL, readID=NULL, 
                      grouping=NULL, psl.rd=NULL, parallel=TRUE) {
  if(!parallel) { registerDoSEQ() }
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(value))
    stopifnot(!is.null(readID))
  } else {
    ## find the otuID by clusters ##
    if(!"clusterTopHit" %in% colnames(psl.rd) | 
       !"clusteredPosition" %in% colnames(psl.rd)) {
          stop("The object supplied in psl.rd parameter does not have 'clusterTopHit' or 'clusteredPosition' column in it. 
          Did you run clusterSites() on it?")
    }
    
    good.rows <- psl.rd$clusterTopHit
    value <- psl.rd$clusteredPosition
    posID <- paste0(space(psl.rd),psl.rd$strand)
    if("qName" %in% colnames(psl.rd)) {
    	readID <- psl.rd$qName
    } else if ("Sequence" %in% colnames(psl.rd)) {
    	readID <- psl.rd$Sequence
    } else {
    	stop("No readID type column found in psl.rd object.")
    }
    
    grouping <- if(is.null(grouping)) { rep("A",nrow(psl.rd)) } else { grouping }
    
    otus <- otuSites2(posID=posID[good.rows], 
                     value=value[good.rows], 
                     readID=readID[good.rows], 
                     grouping=grouping[good.rows], 
                     parallel=parallel)
    
    message("Adding otuIDs back to psl.rd.")        
    otuIDs <- with(otus,
                   split(otuID,
                         paste0(posID,value,readID,grouping)))
    psl.rd$otuID <- NA
    psl.rd$otuID[good.rows] <- as.numeric(otuIDs[paste0(posID, value, readID, 
    													grouping)[good.rows]])
    
    message("Cleaning up!")
    rm("otus","otuIDs","value","posID","readID","grouping","good.rows")
    cleanit <- gc()
    return(psl.rd)
  }
  
  groups <- if(is.null(grouping)) { "A" } else { grouping }
  sites <- data.frame(posID, value, readID,
                      posID2=paste0(posID, value),
                      grouping=groups, stringsAsFactors=FALSE)
  rm(groups)
  
  ## get unique positions per readID by grouping 
  reads <- ddply(sites, .(grouping,readID), summarise, 
                 posIDs=paste(unique(posID2),collapse=","), 
                 counts=length(unique(posID2)), .parallel=parallel)
  
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
  sites <- arrange(sites, grouping, posID, value)
  sites$readID <- as.character(sites$readID)
  sites$grouping <- as.character(sites$grouping)
  
  sites.gr <- with(sites, GRanges(seqnames=posID, IRanges(start=value,width=1), 
                                  strand="*", readID, grouping, counts, 
                                  otuID, newotuID=otuID, check=TRUE))
  mcols(sites.gr)$grouping <- as.character(mcols(sites.gr)$grouping)
  mcols(sites.gr)$readID <- as.character(mcols(sites.gr)$readID)
  
  ## see if readID with a unique/single location matches up to a readID with >1 location, if yes then merge
  mcols(sites.gr)$singles <- mcols(sites.gr)$counts==1
  if(any(mcols(sites.gr)$singles)) {
    message('Merging non-singletons with singletons if any...')    
    sites.gr.list <- split(sites.gr, mcols(sites.gr)$grouping)
    sites.gr <- foreach(x=iter(sites.gr.list), .inorder=FALSE, 
                        .packages="GenomicRanges", .combine=c) %dopar% {
      sigs <- subset(x, mcols(x)$singles)
      nonsigs <- subset(x, !mcols(x)$singles)
      res <- findOverlaps(nonsigs,sigs,maxgap=1)
      if(length(res)>0) {
        res <- as.data.frame(res)
        res$sigsOTU <- mcols(sigs)$otuID[res$subjectHits]
        res$sigsReadID <- mcols(sigs)$readID[res$subjectHits]
        res$nonsigsReadID <- mcols(nonsigs)$readID[res$queryHits]
        s.to.q <- with(res,structure(sigsOTU,names=nonsigsReadID))
        rows <- mcols(nonsigs)$readID %in% res$nonsigsReadID
        mcols(nonsigs)[rows,"newotuID"] <- s.to.q[mcols(nonsigs)[rows,"readID"]]
        mcols(nonsigs)[mcols(nonsigs)$readID %in% res$nonsigsReadID,"check"] <- FALSE
        mcols(sigs)[mcols(sigs)$readID %in% res$sigsReadID,"check"] <- FALSE
      }
      c(sigs,nonsigs)
    }
    rm(sites.gr.list)
  }  
  mcols(sites.gr)$singles <- NULL
  
  ## see if readIDs with >1 locations overlap with other readIDs of the same type ##
  ## this is useful when no readIDs were found with a unique or single locations ##
  ## merge OTUs with overlapping positions within same grouping ##
  message('Merging non-singletons...')
  goods <- subset(sites.gr, !mcols(sites.gr)$check)
  sites.gr <- subset(sites.gr, mcols(sites.gr)$check)
  sites.gr.list <- split(sites.gr, mcols(sites.gr)$grouping)
  sites.gr <- foreach(x=iter(sites.gr.list), .inorder=FALSE, .packages="GenomicRanges", .combine=c) %dopar% {		    
    mcols(x)$readID <- as.character(mcols(x)$readID)
    res <- findOverlaps(x, maxgap=1, ignoreSelf=TRUE,ignoreRedundant=TRUE, select="all")
    if(length(res)>0) {
      res <- as.data.frame(res)
      res$queryOTU <- mcols(x)$otuID[res$queryHits]
      res$queryReadID <- mcols(x)$readID[res$queryHits]       
      res$subjectReadID <- mcols(x)$readID[res$subjectHits]
      s.to.q <- with(res,structure(queryOTU,names=subjectReadID))
      rows <- mcols(x)$readID %in% res$subjectReadID
      mcols(x)[rows,"newotuID"] <- s.to.q[mcols(x)[rows,"readID"]]
      mcols(x)[mcols(x)$readID %in% c(res$subjectReadID,res$queryReadID),"check"] <- FALSE
    }
    x
  }
  sites.gr <- c(sites.gr,goods)
  rm("sites.gr.list","goods")
  cleanit <- gc()
  
  ## trickle the OTU ids back to sites frame ##    
  ots.ids <- sapply(split(mcols(sites.gr)$newotuID,paste0(mcols(sites.gr)$readID, mcols(sites.gr)$grouping)),unique)
  if(!is.numeric(ots.ids)) {
    stop("Something went wrong merging non-singletons. Multiple OTUs assigned to one readID most likely!")
  }
  sites$otuID <- as.numeric(unlist(ots.ids[with(sites,paste0(readID,grouping))],use.names=F))
  
  stopifnot(any(!is.na(sites$otuID)))
  cleanit <- gc()
  
  if(is.null(grouping)) { sites$grouping<-NULL }
  return(sites)
}

#' Remove values/positions which are overlapping between discrete groups based on their frequency.
#'
#' Given a group of discrete factors (i.e. position ids) and integer values, the function tests if they overlap between groups. If overlap is found, then the group having highest frequency of a given position wins, else the position is filtered out from all the groups. The main use of this function is to remove crossover sites from different samples in the data.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand). Required if psl.rd parameter is not defined.
#' @param value a vector of integer locations/positions that needs to be binned, i.e. genomic location. Required if psl.rd parameter is not defined. 
#' @param grouping additional vector of grouping by which to pool the rows (i.e. samplenames). Default is NULL.
#' @param weight a numeric vector of weights to use when calculating frequency of value by posID and grouping if specified. Default is NULL.
#' @param psl.rd a RangedData object. Default is NULL. 
#'
#' @return a data frame of the original input with columns denoting whether a given row is crossover or not. If psl.rd parameter is defined, then a RangedData object with 'isCrossover' column appended at the end.
#'
#' @seealso  \code{\link{clusterSites}}, \code{\link{otuSites}}, \code{\link{otuSites2}}, \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' #crossOverCheck(posID=c('chr1-','chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), value=c(rep(1000,3),5832,1000,12324,65738,928042), grouping=c('a','a','b','b','b','b','c','c'))
#' #crossOverCheck(psl.rd=test.psl.rd)
#'
crossOverCheck <- function(posID=NULL, value=NULL, grouping=NULL, weight=NULL, psl.rd=NULL) {
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(value))
  } else {     
    if(!"clusterTopHit" %in% colnames(psl.rd) | 
         !"clusteredPosition" %in% colnames(psl.rd)) {
      stop("The object supplied in psl.rd parameter does not have 'clusterTopHit' or 'clusteredPosition' column in it. Did you run clusterSites() on it?")
    }
    
    posID <- paste0(space(psl.rd),psl.rd$strand)
    if("clusteredPosition" %in% colnames(psl.rd)) {
      message("Using clusteredPosition column from psl.rd as the value parameter.")
      value <- psl.rd$clusteredPosition
      good.row <- psl.rd$clusterTopHit
    } else if("Position" %in% colnames(psl.rd)) {
      message("Using Position column from psl.rd as the value parameter.")
      value <- psl.rd$Position
      good.row <- rep(TRUE,length(value))
    } else {
      message("Using start(psl.rd) as the value parameter.")
      value <- start(psl.rd)
      good.row <- rep(TRUE,length(value))
    }
    
    isthere <- grepl("isMultiHit",colnames(psl.rd),ignore.case=TRUE)
    if(any(isthere)) { ## see if multihit column exists
      message("Found 'isMultiHit' column in the data. These rows will be ignored for the calculation.")
      good.row <- good.row & !psl.rd[[which(isthere)]]
    }
    
    if(is.null(weight)) { ## see if clonecount column exists
      isthere <- grepl("clonecount",colnames(psl.rd))
      if(any(isthere)) { weight <- psl.rd[[which(isthere)]] } else { weight <- weight }
    }
    grouping <- if(is.null(grouping)) { "" } else { grouping }
    crossed <- crossOverCheck(posID[good.row], value[good.row], grouping=grouping[good.row], weight=weight[good.row])
    
    message("Adding isCrossover data back to psl.rd.")        
    crossedValues <- with(crossed,split(isCrossover,paste0(posID,value,grouping)))
    if(any(sapply(crossedValues,length)>1)) {
    	stop("Error in crossOverCheck: sampling culprits... ", 
    		 paste(names(crossedValues[which(sapply(crossedValues,length)>1)]), collapse=", "))
    }
    psl.rd$isCrossover <- FALSE
    psl.rd$isCrossover[good.row] <- as.logical(crossedValues[paste0(posID,value,grouping)[good.row]])
    
    message("Cleaning up!")
    rm("posID","value","grouping","crossed","crossedValues")
    cleanit <- gc()
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
  res <- findOverlaps(sites.gr, maxgap=1, ignoreSelf=TRUE, ignoreRedundant=FALSE, select="all")
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
#' Given a SampleInfo object, the function finds integration sites for each sample using their respective settings and adds the results back to the object. This is an all-in-one function which BLATs, parses PSL files, finds best hit per read per sample, cluster sites, and assign OTU IDs. It calls \code{\link{blatSeqs}}, \code{\link{read.psl}}, \code{\link{getIntegrationSites}}, \code{\link{clusterSites}}, \code{\link{otuSites}}. There must be linkered reads within the sampleInfo object in order to use this function using the default parameters. If you are planning on BLATing non-linkered reads, then change the seqType to one of accepted options for the 'feature' parameter of \code{\link{extractSeqs}}, except for '!' based features.
#'
#' @param sampleInfo sample information SimpleList object outputted from \code{\link{findLinkers}}, which holds decoded, primed, LTRed, and Linkered sequences for samples per sector/quadrant along with metadata.
#' @param seqType which type of sequence to BLAT and find integration sites for. Default is NULL and determined automatically based on type of restriction enzyme or isolation method used. Could be any one of following options: genomic, genomicLinkered, decoded, primed, LTRed, linkered.
#' @param port a port number to host the gfServer with. Default is 5560.
#' @param host name of the machine running gfServer. Default is 'localhost'.
#' @param genomeIndices an associative character vector of freeze to full or relative path of respective .nib or .2bit files. Default is c("hg18"="/usr/local/blatSuite34/hg18.2bit", "mm8"="/usr/local/blatSuite34/mm8.2bit").
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to TRUE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param samplenames a vector of samplenames to process. Default is NULL, which processes all samples from sampleInfo object.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new data added under each sector and sample. New data attributes include: psl, and sites. The psl attributes holds the BLAT hits per read along with QC information. The sites attribute holds the condensed integration sites where BLAT hits have been clustered by the Position column and cherry picked to have each site pass all the QC steps. 
#'
#' @note If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{findPrimers}}, \code{\link{findLTRs}}, \code{\link{findLinkers}}, \code{\link{startgfServer}}, \code{\link{read.psl}}, \code{\link{blatSeqs}}, \code{\link{blatListedSet}}, \code{\link{pslToRangedObject}}, \code{\link{clusterSites}}, \code{\link{otuSites2}}, \code{\link{crossOverCheck}}, \code{\link{getIntegrationSites}}
#'
#' @export
#'
#' @examples 
#' #findIntegrations(sampleInfo)
#'
findIntegrations <- function(sampleInfo, seqType=NULL, port=5560, host="localhost", genomeIndices=c("hg18"="/usr/local/blatSuite34/hg18.2bit", "mm8"="/usr/local/blatSuite34/mm8.2bit"), parallel=TRUE, samplenames=NULL) {    
    stopifnot(class(sampleInfo)=="SimpleList")
    
    if(!parallel) { registerDoSEQ() }
    
    ## test if there are linkered sequences in the sampleinfo object if specific feature/seqType is not defined ##   
    feature <- ifelse(is.null(seqType),"linkered",seqType)
    
    message("Checking for ",feature," reads.")		
	featured <- extractFeature(sampleInfo,feature=feature)
    samplesfeatured <- sapply(featured,names,simplify=FALSE)
    sectorsfeatured <- names(which(sapply(sapply(featured,length),">",0)))
    rm(featured)
    cleanit <- gc()

    if(length(sectorsfeatured)==0) {
        stop("No ",feature," information found in sampleInfo object provided.")
    }

    ## subset specific samples if defined ##
    samplesToProcess <- unlist(samplesfeatured,use.names=F)
    if(!is.null(samplenames)) {
        samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }

	message("Creating hashes of settings for blatting and processing.")

	## setup settings hashes for blatting the genomic seqs by species & processing hits later ##
    for(setting in c("restrictionenzyme", "freeze", "startwithin", "alignratiothreshold", "genomicpercentidentity", "clustersiteswithin", "keepmultihits")) {
        setter <- extractFeature(sampleInfo,samplename=samplesToProcess,feature=setting)
        names(setter) <- NULL
        setter <- unlist(setter)
        assign(setting,setter)
    }
    
    ## BLAT by respective species ##
    pslFiles <- c()        
    for (f in unique(freeze)) {
        message("Blatting to: ",f)
                
        message("Getting sequences to BLAT")        
        # get sequences to blat #
        if (is.null(seqType)) {
            wanted <- names(restrictionenzyme[!grepl("FRAG|SONIC",restrictionenzyme,ignore.case=TRUE)])
            wanted <- wanted[wanted %in% names(freeze[freeze==f])]
            seqs <- extractSeqs(sampleInfo, samplename=wanted, feature="genomic", minReadLength=5)
            if(any(as.numeric(sapply(seqs,length))>0)) {
                write.listedDNAStringSet(seqs,filePrefix=paste("processed",f,sep=""))
            }

            wanted <- names(restrictionenzyme[grepl("FRAG|SONIC",restrictionenzyme,ignore.case=TRUE)])
            wanted <- wanted[wanted %in% names(freeze[freeze==f])]
            seqs <- extractSeqs(sampleInfo, samplename=wanted, feature="genomicLinkered", minReadLength=5)
            if(any(as.numeric(sapply(seqs,length))>0)) {
                write.listedDNAStringSet(seqs,filePrefix=paste("processed",f,sep=""))
            }
        } else {
            seqs <- extractSeqs(sampleInfo, samplename=names(freeze[freeze==f]), feature=seqType, minReadLength=5)
            if(any(as.numeric(sapply(seqs,length))>0)) {
                write.listedDNAStringSet(seqs,filePrefix=paste("processed",f,sep=""))
            }
        }
        
        # BLAT seqs #
        pslFile <- blatSeqs(query=paste0("processed",f,".*.fa$"), subject=genomeIndices[[f]], standaloneBlat=FALSE, host=host, port=port, parallel=parallel, gzipResults=TRUE)                
        
        message("Cleaning!")
        # add pslFiles for later use #
        pslFiles <- c(pslFiles,pslFile)
        cleanit <- gc()
        system(paste("rm",paste0("processed",f,".*.fa")))
    }
    
    message("Reading PSL files.")
    ## read all hits and split by samples ##
    psl <- read.psl(pslFiles, bestScoring=TRUE, asRangedData=TRUE, removeFile=TRUE, parallel=FALSE)
    cleanit <- gc()
    psl$setname <- sub("^(.+)-(.+)$","\\1",psl$qName)
    psl <- split(psl, psl$setname)
                
    ## begin processing hits ##
    psl.hits <- foreach(x=iter(names(psl)),.inorder=TRUE,.export=c("psl", "sampleInfo", "startwithin", "alignratiothreshold", "genomicpercentidentity", "clustersiteswithin", "keepmultihits")) %dopar% {
        message("Processing ",x)
        
        # add qc info for bestscoring hits #
        psl.x <- getIntegrationSites(psl[[x]], startWithin=startwithin[[x]], alignRatioThreshold=alignratiothreshold[[x]], genomicpercentidentity=genomicpercentidentity[[x]])
        
        # filter multihits if applicable #
        if(!as.logical(keepmultihits[[x]])) {
            psl.x <- psl.x[!psl.x$isMultiHit, ]
        }
        
        # cluster sites by positions #
        psl.x <- clusterSites(psl.rd=psl.x, windowSize=clustersiteswithin[[x]])
        
        # get sites OTU for tagging multihits #
        if(as.logical(keepmultihits[[x]])) {        
            psl.x <- otuSites2(psl.rd=psl.x)
        }
        
        psl.x
    }
    names(psl.hits) <- names(psl)
    
    message("Adding PSL hits back to the object.")
    sampleInfo <- addFeature(sampleInfo,sector=NULL,samplename=names(psl.hits),feature="psl",value=psl.hits)

    message("Adding sites back to the object.")
    psl.hits <- sapply(psl.hits, function(x) {
      x <- DataFrame(x)      
      x <- subset(x, clusterTopHit & pass.allQC, select=setdiff(colnames(x), c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'tSize', 'blockCount', 'blockSizes', 'qStarts', 'score', 'tStarts', 'pass.startWithin', 'alignRatio', 'pass.alignRatio', 'percIdentity', 'pass.percIdentity', 'pass.allQC', 'clusterTopHit', 'width', 'Position')))
      x$start <- x$end <- x$clusteredPosition
      x$clusteredPosition <- NULL    
      RangedData(x)		
    })
    sampleInfo <- addFeature(sampleInfo,sector=NULL,samplename=names(psl.hits),feature="sites",value=psl.hits)

    cleanit <- gc()
        
    sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
    return(sampleInfo)
}