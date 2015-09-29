%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Using hiReadsProcessor}

## Introduction
`hiReadsProcessor` contains set of functions which allow users to process LM-PCR products sequenced using any platform. Given an excel/txt file containing parameters for demultiplexing and sample metadata, the functions in here automate trimming of adaptors, removing any vector sequences, and identifying host genomic location.

The basic philosophy of this package is to detect various 'bits' of information in a sequencing read like barcodes, primers, linkers, etc and store it within a nested list in a space efficient manner. No information is duplicated once read into the object and the supplied utility functions enable addition & extraction of needed features on the fly.

Here is the workflow in a nutshell followed by few simple steps to get you started.

![workflow](hiReadsProcessor_workflow.png)

Please refer to following publications from [Bushman Lab](http://bushmanlab.org/publications) to obtain more information on sample/amplicon preparation

## Quick *hiReadsProcessor* tutorial.

### Step 1: Loading the tools
First load this package and the parallel backend of choice from `BiocParallel` library. Although `BiocParallel` is imported internally and will invoke `multicore` like functionality, you may have to load it to invoke `snow` like functionality.

```r
library(hiReadsProcessor)
## added to avoid package build failure on windows machines ##
if(.Platform$OS.type == "windows") { register(SerialParam()) }
```

### Step 2: Processing the data
The package comes with an example 454 sequencing run: `FLX_sample_run`. In rest of this tutorial we will use this dataset to introduce functionality of this package.

#### Read the data + metadata files

A typical sequencing run will have several filetypes but among the most important files are the fasta/fastq which are created per sector/quadrant/lane. In the example dataset we have compressed fasta files for three quadrants, an excel spreadsheet holding sample metadata along with information/parameters to process that sample and vector fasta files to be trimmed if needed.


```r
runData <- system.file("extdata/FLX_sample_run/", package = "hiReadsProcessor")
list.files(runData, recursive  = TRUE)
```

```
## [1] "RunData/1.TCA.454Reads.fna.gz" "RunData/2.TCA.454Reads.fna.gz"
## [3] "RunData/3.TCA.454Reads.fna.gz" "Vectors/HIV1.fa"              
## [5] "Vectors/MLV-vector.fa"         "sampleInfo.xls"
```

Function `read.SeqFolder` will initiate a `SimpleList` sample information object which will hold everything regarding the sequencing run. The object is structured to store sequencing file paths, sequence data, processed data as well as sample metadata. The function finds the required file needed to ease the automation process. It is important that somewhere in the sequencing folder there is a file called "sampleInfo" else object initialization will fail.


```r
seqProps <- read.SeqFolder(runData, seqfilePattern=".+fna.gz$")
```

```
## Choosing /home/microarray/Bioconductor/hiReadsProcessor/inst/extdata/FLX_sample_run/sampleInfo.xls as sample information file.
```

```r
seqProps
```

```
## List of length 6
## names(6): sequencingFolderPath seqFilePaths ... sectors callHistory
```

On successfully initializing the sample information object you will see that there is a hierarchy within the object. The root or top most level holds the information regarding the sequence folder. Each sample within a quadrant/lane is held within "sectors" list. Within the sectors list, there is a list of samples and the associated metadata which gets modified and appended to as the reads get trimmed and aligned.


```r
seqProps$sectors$"1"$samples
```

```
## List of length 54
## names(54): Roth-MLV-CD4T-20100723NCMu ... Roth-MLV3p-CD4TMLVwell6-MuA
```

#### Demultiplex reads

The first step of processing most sequencing run is to demultiplex reads by barcodes/MIDs. Function `findBarcodes` automates the demultiplexing process based on already stored data from the sample information file. Please see the documentation for `read.sampleInfo` for the kinds of parameters and information held in a sample information file. An example file is supplied within the `FLX_sample_run` dataset.


```r
seqProps <- findBarcodes(seqProps, sector="all", showStats=TRUE)
```

```
## Decoding sector: 1
## Reading:
## /home/microarray/Bioconductor/hiReadsProcessor/inst/extdata/FLX_sample_run/RunData/1.TCA.454Reads.fna.gz
## Using following schema for barcode to sample associations
```

```
##                                barcodesSample
## TGCATCGA           Roth-MLV-CD4T-20100723NCMu
## TCAGTCAG       Roth-MLV-CD4T-20100723well1BMu
## ACGTACGA     Roth-MLV-CD4T-20100723well2BMu1-
## CTCAGACA     Roth-MLV-CD4T-20100723well2BMu10
## CTGAGTCA     Roth-MLV-CD4T-20100723well2BMu11
## TCACAGAC     Roth-MLV-CD4T-20100723well2BMu12
## ACTCGACA      Roth-MLV-CD4T-20100723well2BMu2
## AGACAGTG      Roth-MLV-CD4T-20100723well2BMu3
## CACAGACG      Roth-MLV-CD4T-20100723well2BMu4
## CACTGTGA      Roth-MLV-CD4T-20100723well2BMu5
## CATCTCGA      Roth-MLV-CD4T-20100723well2BMu6
## CATGACGA      Roth-MLV-CD4T-20100723well2BMu7
## CGATCGTA      Roth-MLV-CD4T-20100723well2BMu8
## CGTAGCTA      Roth-MLV-CD4T-20100723well2BMu9
## TCAGTCTC       Roth-MLV-CD4T-20100723well2Mu1
## TCGTAGCA       Roth-MLV-CD4T-20100723well2Mu2
## TCGTCATC       Roth-MLV-CD4T-20100723well2Mu3
## TCTCACAC       Roth-MLV-CD4T-20100723well2Mu4
## TGACAGTC       Roth-MLV-CD4T-20100723well2Mu5
## TGCAGTAC       Roth-MLV-CD4T-20100723well2Mu6
## TCACGTGA           Roth-MLV3p-CD4T-20100730NC
## TGCTGATG   Roth-MLV3p-CD4T-20100730Well1BstYI
## CGTGCGAC    Roth-MLV3p-CD4T-20100730Well1MseI
## CAGCTGTA  Roth-MLV3p-CD4T-20100730Well1NlaIII
## CAGTCTCA Roth-MLV3p-CD4T-20100730Well1Tsp509I
## AGCTCATG   Roth-MLV3p-CD4T-20100730Well2BstyI
## ACACTGAC    Roth-MLV3p-CD4T-20100730Well2MseI
## ACTGAGTC  Roth-MLV3p-CD4T-20100730Well2NlaIII
## TCGATCGA Roth-MLV3p-CD4T-20100730Well2Tsp509I
## AGCACTAC           Roth-MLV3p-CD4TMLVLot16-Mu
## TCTCAGTC         Roth-MLV3p-CD4TMLVWell3-MseI
## TCTCGTCA       Roth-MLV3p-CD4TMLVWell3-NlaIII
## TCGAGTAC      Roth-MLV3p-CD4TMLVWell3-Tsp509I
## TGTGCTGA      Roth-MLV3p-CD4TMLVWell3Harri-Mu
## ACACACTG      Roth-MLV3p-CD4TMLVWell3Lot60-Mu
## ACAGTGTC      Roth-MLV3p-CD4TMLVWell3Lot62-Mu
## ACGATGCT      Roth-MLV3p-CD4TMLVWell3Lot64-Mu
## ACGTCATG   Roth-MLV3p-CD4TMLVWell3Lot64new-Mu
## TGACTGTG        Roth-MLV3p-CD4TMLVWell4-BstYI
## GCTATACA         Roth-MLV3p-CD4TMLVWell4-MseI
## GAGCATGA       Roth-MLV3p-CD4TMLVWell4-NlaIII
## GCATGCTA      Roth-MLV3p-CD4TMLVWell4-Tsp509I
## ACTGATAC        Roth-MLV3p-CD4TMLVWell5-BstYI
## GTCTGAGC         Roth-MLV3p-CD4TMLVWell5-MseI
## AGCTCTGC       Roth-MLV3p-CD4TMLVWell5-NlaIII
## ATACTCTC      Roth-MLV3p-CD4TMLVWell5-Tsp509I
## ACTGACAC        Roth-MLV3p-CD4TMLVWell6-BstYI
## TGACGTCA         Roth-MLV3p-CD4TMLVWell6-MseI
## CAGTCACG       Roth-MLV3p-CD4TMLVWell6-NlaIII
## TCGAGCAT      Roth-MLV3p-CD4TMLVWell6-Tsp509I
## AGCTGTAC        Roth-MLV3p-CD4TMLVwell3-BstYI
## TCGAGACT          Roth-MLV3p-CD4TMLVwell4-MuA
## TCGACTGA          Roth-MLV3p-CD4TMLVwell5-MuA
## ACAGCAGA          Roth-MLV3p-CD4TMLVwell6-MuA
```

```
## Number of Sequences with no matching barcode: 196
## Number of Sequences decoded:
```

```
##                             sampleNames Freq
## 1        Roth-MLV-CD4T-20100723well1BMu  219
## 2      Roth-MLV-CD4T-20100723well2BMu1-   37
## 3      Roth-MLV-CD4T-20100723well2BMu10   12
## 4      Roth-MLV-CD4T-20100723well2BMu11   40
## 5      Roth-MLV-CD4T-20100723well2BMu12    9
## 6       Roth-MLV-CD4T-20100723well2BMu2   50
## 7       Roth-MLV-CD4T-20100723well2BMu3   88
## 8       Roth-MLV-CD4T-20100723well2BMu4   85
## 9       Roth-MLV-CD4T-20100723well2BMu5   86
## 10      Roth-MLV-CD4T-20100723well2BMu6   85
## 11      Roth-MLV-CD4T-20100723well2BMu7   47
## 12      Roth-MLV-CD4T-20100723well2BMu8   26
## 13      Roth-MLV-CD4T-20100723well2BMu9   28
## 14       Roth-MLV-CD4T-20100723well2Mu1   57
## 15       Roth-MLV-CD4T-20100723well2Mu2   65
## 16       Roth-MLV-CD4T-20100723well2Mu3   31
## 17       Roth-MLV-CD4T-20100723well2Mu4   22
## 18       Roth-MLV-CD4T-20100723well2Mu5    6
## 19       Roth-MLV-CD4T-20100723well2Mu6    2
## 20   Roth-MLV3p-CD4T-20100730Well1BstYI  139
## 21    Roth-MLV3p-CD4T-20100730Well1MseI  240
## 22  Roth-MLV3p-CD4T-20100730Well1NlaIII  251
## 23 Roth-MLV3p-CD4T-20100730Well1Tsp509I  313
## 24   Roth-MLV3p-CD4T-20100730Well2BstyI  324
## 25    Roth-MLV3p-CD4T-20100730Well2MseI  475
## 26  Roth-MLV3p-CD4T-20100730Well2NlaIII  371
## 27 Roth-MLV3p-CD4T-20100730Well2Tsp509I  437
## 28           Roth-MLV3p-CD4TMLVLot16-Mu  172
## 29         Roth-MLV3p-CD4TMLVWell3-MseI  468
## 30       Roth-MLV3p-CD4TMLVWell3-NlaIII  301
## 31      Roth-MLV3p-CD4TMLVWell3-Tsp509I  487
## 32      Roth-MLV3p-CD4TMLVWell3Harri-Mu   62
## 33      Roth-MLV3p-CD4TMLVWell3Lot60-Mu   36
## 34      Roth-MLV3p-CD4TMLVWell3Lot62-Mu  171
## 35      Roth-MLV3p-CD4TMLVWell3Lot64-Mu  109
## 36   Roth-MLV3p-CD4TMLVWell3Lot64new-Mu   14
## 37        Roth-MLV3p-CD4TMLVWell4-BstYI  116
## 38         Roth-MLV3p-CD4TMLVWell4-MseI  395
## 39       Roth-MLV3p-CD4TMLVWell4-NlaIII  287
## 40      Roth-MLV3p-CD4TMLVWell4-Tsp509I  318
## 41        Roth-MLV3p-CD4TMLVWell5-BstYI  148
## 42         Roth-MLV3p-CD4TMLVWell5-MseI  284
## 43       Roth-MLV3p-CD4TMLVWell5-NlaIII  286
## 44      Roth-MLV3p-CD4TMLVWell5-Tsp509I  542
## 45        Roth-MLV3p-CD4TMLVWell6-BstYI  195
## 46         Roth-MLV3p-CD4TMLVWell6-MseI  378
## 47       Roth-MLV3p-CD4TMLVWell6-NlaIII  261
## 48      Roth-MLV3p-CD4TMLVWell6-Tsp509I  490
## 49        Roth-MLV3p-CD4TMLVwell3-BstYI  211
## 50          Roth-MLV3p-CD4TMLVwell4-MuA  138
## 51          Roth-MLV3p-CD4TMLVwell5-MuA  235
## 52          Roth-MLV3p-CD4TMLVwell6-MuA  155
```

```
## Decoding sector: 2
## Reading:
## /home/microarray/Bioconductor/hiReadsProcessor/inst/extdata/FLX_sample_run/RunData/2.TCA.454Reads.fna.gz
## Using following schema for barcode to sample associations
```

```
##                                barcodesSample
## TGCATCGA           Roth-MLV-CD4T-20100723NCMu
## TCAGTCAG       Roth-MLV-CD4T-20100723well1BMu
## ACGTACGA     Roth-MLV-CD4T-20100723well2BMu1-
## CTCAGACA     Roth-MLV-CD4T-20100723well2BMu10
## CTGAGTCA     Roth-MLV-CD4T-20100723well2BMu11
## TCACAGAC     Roth-MLV-CD4T-20100723well2BMu12
## ACTCGACA      Roth-MLV-CD4T-20100723well2BMu2
## AGACAGTG      Roth-MLV-CD4T-20100723well2BMu3
## CACAGACG      Roth-MLV-CD4T-20100723well2BMu4
## CACTGTGA      Roth-MLV-CD4T-20100723well2BMu5
## CATCTCGA      Roth-MLV-CD4T-20100723well2BMu6
## CATGACGA      Roth-MLV-CD4T-20100723well2BMu7
## CGATCGTA      Roth-MLV-CD4T-20100723well2BMu8
## CGTAGCTA      Roth-MLV-CD4T-20100723well2BMu9
## TCAGTCTC       Roth-MLV-CD4T-20100723well2Mu1
## TCGTAGCA       Roth-MLV-CD4T-20100723well2Mu2
## TCGTCATC       Roth-MLV-CD4T-20100723well2Mu3
## TCTCACAC       Roth-MLV-CD4T-20100723well2Mu4
## TGACAGTC       Roth-MLV-CD4T-20100723well2Mu5
## TGCAGTAC       Roth-MLV-CD4T-20100723well2Mu6
## TCACGTGA           Roth-MLV3p-CD4T-20100730NC
## TGCTGATG   Roth-MLV3p-CD4T-20100730Well1BstYI
## CGTGCGAC    Roth-MLV3p-CD4T-20100730Well1MseI
## CAGCTGTA  Roth-MLV3p-CD4T-20100730Well1NlaIII
## CAGTCTCA Roth-MLV3p-CD4T-20100730Well1Tsp509I
## AGCTCATG   Roth-MLV3p-CD4T-20100730Well2BstyI
## ACACTGAC    Roth-MLV3p-CD4T-20100730Well2MseI
## ACTGAGTC  Roth-MLV3p-CD4T-20100730Well2NlaIII
## TCGATCGA Roth-MLV3p-CD4T-20100730Well2Tsp509I
## AGCACTAC           Roth-MLV3p-CD4TMLVLot16-Mu
## TCTCAGTC         Roth-MLV3p-CD4TMLVWell3-MseI
## TCTCGTCA       Roth-MLV3p-CD4TMLVWell3-NlaIII
## TCGAGTAC      Roth-MLV3p-CD4TMLVWell3-Tsp509I
## TGTGCTGA      Roth-MLV3p-CD4TMLVWell3Harri-Mu
## ACACACTG      Roth-MLV3p-CD4TMLVWell3Lot60-Mu
## ACAGTGTC      Roth-MLV3p-CD4TMLVWell3Lot62-Mu
## ACGATGCT      Roth-MLV3p-CD4TMLVWell3Lot64-Mu
## ACGTCATG   Roth-MLV3p-CD4TMLVWell3Lot64new-Mu
## TGACTGTG        Roth-MLV3p-CD4TMLVWell4-BstYI
## GCTATACA         Roth-MLV3p-CD4TMLVWell4-MseI
## GAGCATGA       Roth-MLV3p-CD4TMLVWell4-NlaIII
## GCATGCTA      Roth-MLV3p-CD4TMLVWell4-Tsp509I
## ACTGATAC        Roth-MLV3p-CD4TMLVWell5-BstYI
## GTCTGAGC         Roth-MLV3p-CD4TMLVWell5-MseI
## AGCTCTGC       Roth-MLV3p-CD4TMLVWell5-NlaIII
## ATACTCTC      Roth-MLV3p-CD4TMLVWell5-Tsp509I
## ACTGACAC        Roth-MLV3p-CD4TMLVWell6-BstYI
## TGACGTCA         Roth-MLV3p-CD4TMLVWell6-MseI
## CAGTCACG       Roth-MLV3p-CD4TMLVWell6-NlaIII
## TCGAGCAT      Roth-MLV3p-CD4TMLVWell6-Tsp509I
## AGCTGTAC        Roth-MLV3p-CD4TMLVwell3-BstYI
## TCGAGACT          Roth-MLV3p-CD4TMLVwell4-MuA
## TCGACTGA          Roth-MLV3p-CD4TMLVwell5-MuA
## ACAGCAGA          Roth-MLV3p-CD4TMLVwell6-MuA
```

```
## Number of Sequences with no matching barcode: 214
## Number of Sequences decoded:
```

```
##                             sampleNames Freq
## 1        Roth-MLV-CD4T-20100723well1BMu  215
## 2      Roth-MLV-CD4T-20100723well2BMu1-   27
## 3      Roth-MLV-CD4T-20100723well2BMu10    4
## 4      Roth-MLV-CD4T-20100723well2BMu11   41
## 5      Roth-MLV-CD4T-20100723well2BMu12    7
## 6       Roth-MLV-CD4T-20100723well2BMu2   41
## 7       Roth-MLV-CD4T-20100723well2BMu3   78
## 8       Roth-MLV-CD4T-20100723well2BMu4   76
## 9       Roth-MLV-CD4T-20100723well2BMu5   74
## 10      Roth-MLV-CD4T-20100723well2BMu6   80
## 11      Roth-MLV-CD4T-20100723well2BMu7   47
## 12      Roth-MLV-CD4T-20100723well2BMu8   15
## 13      Roth-MLV-CD4T-20100723well2BMu9   32
## 14       Roth-MLV-CD4T-20100723well2Mu1   50
## 15       Roth-MLV-CD4T-20100723well2Mu2   68
## 16       Roth-MLV-CD4T-20100723well2Mu3   21
## 17       Roth-MLV-CD4T-20100723well2Mu4   19
## 18       Roth-MLV-CD4T-20100723well2Mu5    4
## 19       Roth-MLV-CD4T-20100723well2Mu6    2
## 20   Roth-MLV3p-CD4T-20100730Well1BstYI  151
## 21    Roth-MLV3p-CD4T-20100730Well1MseI  223
## 22  Roth-MLV3p-CD4T-20100730Well1NlaIII  253
## 23 Roth-MLV3p-CD4T-20100730Well1Tsp509I  315
## 24   Roth-MLV3p-CD4T-20100730Well2BstyI  375
## 25    Roth-MLV3p-CD4T-20100730Well2MseI  428
## 26  Roth-MLV3p-CD4T-20100730Well2NlaIII  370
## 27 Roth-MLV3p-CD4T-20100730Well2Tsp509I  456
## 28           Roth-MLV3p-CD4TMLVLot16-Mu  153
## 29         Roth-MLV3p-CD4TMLVWell3-MseI  487
## 30       Roth-MLV3p-CD4TMLVWell3-NlaIII  311
## 31      Roth-MLV3p-CD4TMLVWell3-Tsp509I  469
## 32      Roth-MLV3p-CD4TMLVWell3Harri-Mu   54
## 33      Roth-MLV3p-CD4TMLVWell3Lot60-Mu   54
## 34      Roth-MLV3p-CD4TMLVWell3Lot62-Mu  165
## 35      Roth-MLV3p-CD4TMLVWell3Lot64-Mu  116
## 36   Roth-MLV3p-CD4TMLVWell3Lot64new-Mu    9
## 37        Roth-MLV3p-CD4TMLVWell4-BstYI  111
## 38         Roth-MLV3p-CD4TMLVWell4-MseI  400
## 39       Roth-MLV3p-CD4TMLVWell4-NlaIII  286
## 40      Roth-MLV3p-CD4TMLVWell4-Tsp509I  338
## 41        Roth-MLV3p-CD4TMLVWell5-BstYI  170
## 42         Roth-MLV3p-CD4TMLVWell5-MseI  319
## 43       Roth-MLV3p-CD4TMLVWell5-NlaIII  283
## 44      Roth-MLV3p-CD4TMLVWell5-Tsp509I  608
## 45        Roth-MLV3p-CD4TMLVWell6-BstYI  205
## 46         Roth-MLV3p-CD4TMLVWell6-MseI  406
## 47       Roth-MLV3p-CD4TMLVWell6-NlaIII  247
## 48      Roth-MLV3p-CD4TMLVWell6-Tsp509I  445
## 49        Roth-MLV3p-CD4TMLVwell3-BstYI  179
## 50          Roth-MLV3p-CD4TMLVwell4-MuA  125
## 51          Roth-MLV3p-CD4TMLVwell5-MuA  250
## 52          Roth-MLV3p-CD4TMLVwell6-MuA  124
```

```
## Decoding sector: 3
## Reading:
## /home/microarray/Bioconductor/hiReadsProcessor/inst/extdata/FLX_sample_run/RunData/3.TCA.454Reads.fna.gz
## Using following schema for barcode to sample associations
```

```
##                                 barcodesSample
## CCGGAATT   Ocwieja-HIV896-CD4TND365-InfectionI
## GATCGACT  Ocwieja-HIV896-CD4TND365-InfectionII
## TCGTACAG Ocwieja-HIV896-CD4TND365-InfectionIII
## TATAGCGC      Ocwieja-HIV896-CD4TND365-NoVirus
```

```
## Number of Sequences with no matching barcode: 140
## Number of Sequences decoded:
```

```
##                             sampleNames Freq
## 1   Ocwieja-HIV896-CD4TND365-InfectionI 3879
## 2  Ocwieja-HIV896-CD4TND365-InfectionII 3111
## 3 Ocwieja-HIV896-CD4TND365-InfectionIII 2869
## 4      Ocwieja-HIV896-CD4TND365-NoVirus    1
```

```r
seqProps
```

```
## List of length 6
## names(6): sequencingFolderPath seqFilePaths ... sectors callHistory
```

```r
seqProps$sectors
```

```
## List of length 3
## names(3): 1 2 3
```

#### Detect primers

Following the barcode sequence is the 5' viral LTR primer. Function `findPrimers` facilities the trimming of respective primers(__primerltrsequence__) for each sample. Minimum threshold for detecting the primer can be adjusted using __primerLTRidentity__ within the sample information file.

Since it may take a while to process this kind of data, the package is equipped with processed data object which makes things easier for the tutorial. The code chunks below are not evaluated but rather references the loaded `seqProps` object.


```r
load(file.path(system.file("data", package = "hiReadsProcessor"),
               "FLX_seqProps.RData"))
```


```r
seqProps <- findPrimers(seqProps, showStats=TRUE)
```

#### Detect LTR edge

If LM-PCR products were designed to include the viral LTR (__ltrBitSequence__) following the primer landing site, then `findLTRs` confirms authenticity of the integrated virus. Absence of LTR part denotes nongenuine integration! Minimum threshold for detecting the LTR bit can be adjusted using __ltrBitIdentity__.


```r
seqProps <- findLTRs(seqProps, showStats=TRUE)
```

#### Detect vector sequence...if any!

If the __vectorFile__ parameter is defined within the sample information file, function `findVector` tags any reads which matches the given vector file. These reads are discarded during the genomic alignment step which is covered later. 


```r
seqProps <- findVector(seqProps, showStats=TRUE)
```

#### Detect Linker adaptors

Linker adaptors are found on the 3' end of sequences. Depending on an experiment the __linkerSequence__ can be same or different per sample. Furthermore, some linker adaptors are designed to have __*primerID*__ which can help quantify pre-PCR products. Function `findLinkers` makes it easy to process various samples with different linker sequences and type. If primerID technology is utilized, enabling parameter __primerIdInLinker__ within the sample information file automates the extract of the random part within the adaptor. Thresholds for linker detection can be controlled by setting following parameters within the sample information file: __linkerIdentity, primerIdInLinker, primerIdInLinkerIdentity1, primerIdInLinkerIdentity2__


```r
seqProps <- findLinkers(seqProps, showStats=TRUE, doRC=TRUE)
```

#### Detect integration sites

Once all the non-genomic parts have been detected, it is time to find the actual integration sites. Function `findIntegrations` makes this a breeze given that BLAT and indexed genome files are provided/in-place.


```r
seqProps <- findIntegrations(seqProps, 
                             genomeIndices=c("hg18"="/usr/local/genomeIndexes/hg18.noRandom.2bit"), numServers=2)
```

#### Obtain summary

Function `sampleSummary` quantifies 7 basic features of this package: 

* decoded: # of reads which found an exact match to the barcode
* primed: # of reads with primer detected after the barcode
* LTRed: # of reads with LTR bit detected after the primer
* vectored: # of reads with vector sequence following the LTR bit
* linkered: # of reads with linker detected
* psl: # of reads with primer, LTR bit, and/or linker which aligned to the genome
* sites: # of unique integration sites (dereplicated by genomic position)


```r
sampleSummary(seqProps)
```

```
## Total sectors:1,2,3
```

```
## Warning in rbind_all(x, .id): Unequal factor levels: coercing to character
```

```
## Source: local data frame [112 x 9]
## 
##    Sector                       SampleName decoded primed LTRed vectored
##     (chr)                            (chr)   (int)  (int) (int)    (int)
## 1       1       Roth-MLV-CD4T-20100723NCMu      NA     NA    NA       NA
## 2       1   Roth-MLV-CD4T-20100723well1BMu      18     18     6        1
## 3       1 Roth-MLV-CD4T-20100723well2BMu1-      12     11    NA       NA
## 4       1 Roth-MLV-CD4T-20100723well2BMu10       1      1     1       NA
## 5       1 Roth-MLV-CD4T-20100723well2BMu11       6      4     3       NA
## 6       1 Roth-MLV-CD4T-20100723well2BMu12       1      1    NA       NA
## 7       1  Roth-MLV-CD4T-20100723well2BMu2       3      3     1       NA
## 8       1  Roth-MLV-CD4T-20100723well2BMu3       8      8    NA       NA
## 9       1  Roth-MLV-CD4T-20100723well2BMu4      10     10     4       NA
## 10      1  Roth-MLV-CD4T-20100723well2BMu5       8      8     3       NA
## ..    ...                              ...     ...    ...   ...      ...
## Variables not shown: linkered (int), psl (int), sites (int)
```

## Detailed *hiReadsProcessor* tutorial.

Before diving into functions offered by this package, lets first understand the underlying data object holding all the data. For example purposes we will refer to this as the "sampleInfo" object (although it's essentially a SimpleList object).

The figure below outlines the hierarchy of data storage within the sampleInfo object.
![sampleInfoObj](hiReadsProcessor_object.png)

** THE SECTIONS BELOW ARE IN WORKS **
### Sequencing Run related functions
* read.SeqFolder
* read.sampleInfo
* read.seqsFromSector
* findBarcodes|decodeByBarcode
* findPrimers
* findLTRs
* findVector
* findLinkers
* troubleshootLinkers

### Reads related functions
* splitByBarcode
* dereplicateReads
* replicateReads
* splitSeqsToFiles
* write.listedDNAStringSet
* trimSeqs
* removeReadsWithNs

### Alignment functions
* startgfServer
* stopgfServer
* blatSeqs
* blatListedSet
* subreadAlignSeqs
* vpairwiseAlignSeqs
* pairwiseAlignSeqs
* primerIDAlignSeqs

### Alignment related functions
* read.psl
* write.psl
* read.blast8
* read.BAMasPSL
* pairUpAlignments

### Integration Site functions
* findIntegrations
* getIntegrationSites
* clusterSites
* getSonicAbund
* isuSites, otuSites
* crossOverCheck
* annotateSites

### Utility functions
* findAndTrimSeq
* sampleSummary
* pslToRangedObject
* pslCols
* getSectorsForSamples
* extractSeqs
* extractFeature
* doRCtest
* chunkize
* addFeature
* addListNameToReads

## Session Info

```r
sessionInfo()
```

```
## R version 3.2.2 (2015-08-14)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## Running under: CentOS release 6.7 (Final)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8          LC_NUMERIC=C                 
##  [3] LC_TIME=en_US.UTF-8           LC_COLLATE=C                 
##  [5] LC_MONETARY=en_US.UTF-8       LC_MESSAGES=en_US.UTF-8      
##  [7] LC_PAPER=en_US.UTF-8          LC_NAME=en_US.UTF-8          
##  [9] LC_ADDRESS=en_US.UTF-8        LC_TELEPHONE=en_US.UTF-8     
## [11] LC_MEASUREMENT=en_US.UTF-8    LC_IDENTIFICATION=en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] hiReadsProcessor_1.1.4  hiAnnotator_1.3.2      
##  [3] BiocParallel_1.2.21     xlsx_0.5.7             
##  [5] xlsxjars_0.6.1          rJava_0.9-7            
##  [7] GenomicAlignments_1.4.1 Rsamtools_1.20.4       
##  [9] GenomicRanges_1.20.8    GenomeInfoDb_1.4.2     
## [11] Biostrings_2.36.4       XVector_0.8.0          
## [13] IRanges_2.2.7           S4Vectors_0.6.6        
## [15] BiocGenerics_0.14.0    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.1          formatR_1.2.1        RColorBrewer_1.1-2  
##  [4] futile.logger_1.4.1  plyr_1.8.3           bitops_1.0-6        
##  [7] futile.options_1.0.0 iterators_1.0.7      tools_3.2.2         
## [10] zlibbioc_1.14.0      digest_0.6.8         evaluate_0.8        
## [13] lattice_0.20-33      memoise_0.2.1        gtable_0.1.2        
## [16] BSgenome_1.36.3      foreach_1.4.2        DBI_0.3.1           
## [19] proto_0.3-10         knitr_1.11           hwriter_1.3.2       
## [22] rSFFreader_0.16.0    dplyr_0.4.3          rtracklayer_1.28.10 
## [25] stringr_1.0.0        roxygen2_4.1.1       devtools_1.9.1      
## [28] grid_3.2.2           Biobase_2.28.0       R6_2.1.1            
## [31] XML_3.98-1.3         latticeExtra_0.6-26  ggplot2_1.0.1       
## [34] reshape2_1.4.1       lambda.r_1.1.7       magrittr_1.5        
## [37] splines_3.2.2        MASS_7.3-44          scales_0.3.0        
## [40] codetools_0.2-14     ShortRead_1.26.0     assertthat_0.1      
## [43] colorspace_1.2-6     sonicLength_1.4.4    stringi_0.5-5       
## [46] munsell_0.4.2        RCurl_1.95-4.7
```

