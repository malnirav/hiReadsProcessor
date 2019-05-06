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
##  [1] "RunData/1.TCA.454Reads.fna"    "RunData/1.TCA.454Reads.fna.gz"
##  [3] "RunData/2.TCA.454Reads.fna"    "RunData/2.TCA.454Reads.fna.gz"
##  [5] "RunData/3.TCA.454Reads.fna"    "RunData/3.TCA.454Reads.fna.gz"
##  [7] "Vectors/HIV1.fa"               "Vectors/MLV-vector.fa"        
##  [9] "sampleInfo.xls"                "sampleInfo.xlsx"
```

Function `read.SeqFolder` will initiate a `SimpleList` sample information object which will hold everything regarding the sequencing run. The object is structured to store sequencing file paths, sequence data, processed data as well as sample metadata. The function finds the required file needed to ease the automation process. It is important that somewhere in the sequencing folder there is a file called "sampleInfo" else object initialization will fail.


```r
seqProps <- read.SeqFolder(runData, seqfilePattern = ".+fna.gz$")
```

```
## Choosing /Users/nerv/github/hiReadsProcessor/inst/extdata/FLX_sample_run/sampleInfo.xls as sample information file.
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
seqProps <- findBarcodes(seqProps, sector = "all", showStats = TRUE)
```

```
## Decoding sector: 1
```

```
## Reading:
## /Users/nerv/github/hiReadsProcessor/inst/extdata/FLX_sample_run/RunData/1.TCA.454Reads.fna.gz
```

```
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
```

```
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
```

```
## Reading:
## /Users/nerv/github/hiReadsProcessor/inst/extdata/FLX_sample_run/RunData/2.TCA.454Reads.fna.gz
```

```
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
```

```
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
```

```
## Reading:
## /Users/nerv/github/hiReadsProcessor/inst/extdata/FLX_sample_run/RunData/3.TCA.454Reads.fna.gz
```

```
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
```

```
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
seqProps <- findVector(seqProps, showStats = TRUE)
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
                             genomeIndices = c("hg18" = "/usr/local/genomeIndexes/hg18.noRandom.2bit"), 
                             numServers = 2)
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
## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character
```

```
## Warning in bind_rows_(x, .id): binding character and factor vector,
## coercing into character vector

## Warning in bind_rows_(x, .id): binding character and factor vector,
## coercing into character vector

## Warning in bind_rows_(x, .id): binding character and factor vector,
## coercing into character vector
```

```
##     Sector                            SampleName decoded primed LTRed
## 1        1            Roth-MLV-CD4T-20100723NCMu      NA     NA    NA
## 2        1        Roth-MLV-CD4T-20100723well1BMu      18     18     6
## 3        1      Roth-MLV-CD4T-20100723well2BMu1-      12     11    NA
## 4        1      Roth-MLV-CD4T-20100723well2BMu10       1      1     1
## 5        1      Roth-MLV-CD4T-20100723well2BMu11       6      4     3
## 6        1      Roth-MLV-CD4T-20100723well2BMu12       1      1    NA
## 7        1       Roth-MLV-CD4T-20100723well2BMu2       3      3     1
## 8        1       Roth-MLV-CD4T-20100723well2BMu3       8      8    NA
## 9        1       Roth-MLV-CD4T-20100723well2BMu4      10     10     4
## 10       1       Roth-MLV-CD4T-20100723well2BMu5       8      8     3
## 11       1       Roth-MLV-CD4T-20100723well2BMu6       9      8     2
## 12       1       Roth-MLV-CD4T-20100723well2BMu7       7      5     1
## 13       1       Roth-MLV-CD4T-20100723well2BMu8       4      4     2
## 14       1       Roth-MLV-CD4T-20100723well2BMu9       1      1     1
## 15       1        Roth-MLV-CD4T-20100723well2Mu1       3      3    NA
## 16       1        Roth-MLV-CD4T-20100723well2Mu2       4      4    NA
## 17       1        Roth-MLV-CD4T-20100723well2Mu3       3      3    NA
## 18       1        Roth-MLV-CD4T-20100723well2Mu4       1      1    NA
## 19       1        Roth-MLV-CD4T-20100723well2Mu5       1      1    NA
## 20       1        Roth-MLV-CD4T-20100723well2Mu6      NA     NA    NA
## 21       1            Roth-MLV3p-CD4T-20100730NC      NA     NA    NA
## 22       1    Roth-MLV3p-CD4T-20100730Well1BstYI      19     18    15
## 23       1     Roth-MLV3p-CD4T-20100730Well1MseI      26     26    23
## 24       1   Roth-MLV3p-CD4T-20100730Well1NlaIII      36     33    28
## 25       1  Roth-MLV3p-CD4T-20100730Well1Tsp509I      34     29    25
## 26       1    Roth-MLV3p-CD4T-20100730Well2BstyI      32     31    22
## 27       1     Roth-MLV3p-CD4T-20100730Well2MseI      43     43    39
## 28       1   Roth-MLV3p-CD4T-20100730Well2NlaIII      31     31    25
## 29       1  Roth-MLV3p-CD4T-20100730Well2Tsp509I      39     39    35
## 30       1            Roth-MLV3p-CD4TMLVLot16-Mu      15     14    11
## 31       1         Roth-MLV3p-CD4TMLVwell3-BstYI      18     18    17
## 32       1          Roth-MLV3p-CD4TMLVWell3-MseI      42     40    36
## 33       1        Roth-MLV3p-CD4TMLVWell3-NlaIII      35     35    28
## 34       1       Roth-MLV3p-CD4TMLVWell3-Tsp509I      39     39    37
## 35       1       Roth-MLV3p-CD4TMLVWell3Harri-Mu       9      9     7
## 36       1       Roth-MLV3p-CD4TMLVWell3Lot60-Mu       5      5     2
## 37       1       Roth-MLV3p-CD4TMLVWell3Lot62-Mu      13     13     8
## 38       1       Roth-MLV3p-CD4TMLVWell3Lot64-Mu       8      8     8
## 39       1    Roth-MLV3p-CD4TMLVWell3Lot64new-Mu       2      2     2
## 40       1         Roth-MLV3p-CD4TMLVWell4-BstYI       8      8     6
## 41       1          Roth-MLV3p-CD4TMLVWell4-MseI      38     38    35
## 42       1           Roth-MLV3p-CD4TMLVwell4-MuA      18     18     1
## 43       1        Roth-MLV3p-CD4TMLVWell4-NlaIII      22     22    22
## 44       1       Roth-MLV3p-CD4TMLVWell4-Tsp509I      30     26    25
## 45       1         Roth-MLV3p-CD4TMLVWell5-BstYI      10     10     8
## 46       1          Roth-MLV3p-CD4TMLVWell5-MseI      43     41    38
## 47       1           Roth-MLV3p-CD4TMLVwell5-MuA      26     25     6
## 48       1        Roth-MLV3p-CD4TMLVWell5-NlaIII      28     27    27
## 49       1       Roth-MLV3p-CD4TMLVWell5-Tsp509I      67     66    60
## 50       1         Roth-MLV3p-CD4TMLVWell6-BstYI      17     17    14
## 51       1          Roth-MLV3p-CD4TMLVWell6-MseI      36     36    31
## 52       1           Roth-MLV3p-CD4TMLVwell6-MuA      13     13    NA
## 53       1        Roth-MLV3p-CD4TMLVWell6-NlaIII      22     22    21
## 54       1       Roth-MLV3p-CD4TMLVWell6-Tsp509I      54     52    48
## 55       2            Roth-MLV-CD4T-20100723NCMu      NA     NA    NA
## 56       2        Roth-MLV-CD4T-20100723well1BMu      21     21     6
## 57       2      Roth-MLV-CD4T-20100723well2BMu1-       2      2    NA
## 58       2      Roth-MLV-CD4T-20100723well2BMu10      NA     NA    NA
## 59       2      Roth-MLV-CD4T-20100723well2BMu11       4      4     2
## 60       2      Roth-MLV-CD4T-20100723well2BMu12       1      1    NA
## 61       2       Roth-MLV-CD4T-20100723well2BMu2       4      4    NA
## 62       2       Roth-MLV-CD4T-20100723well2BMu3      17     17     2
## 63       2       Roth-MLV-CD4T-20100723well2BMu4       8      8     1
## 64       2       Roth-MLV-CD4T-20100723well2BMu5       6      6     3
## 65       2       Roth-MLV-CD4T-20100723well2BMu6       5      5     1
## 66       2       Roth-MLV-CD4T-20100723well2BMu7       8      6     2
## 67       2       Roth-MLV-CD4T-20100723well2BMu8       2      2     1
## 68       2       Roth-MLV-CD4T-20100723well2BMu9       3      3     3
## 69       2        Roth-MLV-CD4T-20100723well2Mu1       2      2    NA
## 70       2        Roth-MLV-CD4T-20100723well2Mu2       7      7    NA
## 71       2        Roth-MLV-CD4T-20100723well2Mu3       4      4    NA
## 72       2        Roth-MLV-CD4T-20100723well2Mu4       1      1    NA
## 73       2        Roth-MLV-CD4T-20100723well2Mu5      NA     NA    NA
## 74       2        Roth-MLV-CD4T-20100723well2Mu6      NA     NA    NA
## 75       2            Roth-MLV3p-CD4T-20100730NC      NA     NA    NA
## 76       2    Roth-MLV3p-CD4T-20100730Well1BstYI      14     13     8
## 77       2     Roth-MLV3p-CD4T-20100730Well1MseI      20     20    20
## 78       2   Roth-MLV3p-CD4T-20100730Well1NlaIII      29     28    24
## 79       2  Roth-MLV3p-CD4T-20100730Well1Tsp509I      25     22    20
## 80       2    Roth-MLV3p-CD4T-20100730Well2BstyI      45     45    36
## 81       2     Roth-MLV3p-CD4T-20100730Well2MseI      46     46    45
## 82       2   Roth-MLV3p-CD4T-20100730Well2NlaIII      36     35    30
## 83       2  Roth-MLV3p-CD4T-20100730Well2Tsp509I      34     33    31
## 84       2            Roth-MLV3p-CD4TMLVLot16-Mu      16     16     9
## 85       2         Roth-MLV3p-CD4TMLVwell3-BstYI      17     17    14
## 86       2          Roth-MLV3p-CD4TMLVWell3-MseI      64     63    57
## 87       2        Roth-MLV3p-CD4TMLVWell3-NlaIII      32     32    28
## 88       2       Roth-MLV3p-CD4TMLVWell3-Tsp509I      42     42    37
## 89       2       Roth-MLV3p-CD4TMLVWell3Harri-Mu       3      3     1
## 90       2       Roth-MLV3p-CD4TMLVWell3Lot60-Mu       3      3     3
## 91       2       Roth-MLV3p-CD4TMLVWell3Lot62-Mu      16     16    10
## 92       2       Roth-MLV3p-CD4TMLVWell3Lot64-Mu      11     11    10
## 93       2    Roth-MLV3p-CD4TMLVWell3Lot64new-Mu      NA     NA    NA
## 94       2         Roth-MLV3p-CD4TMLVWell4-BstYI      13     13    13
## 95       2          Roth-MLV3p-CD4TMLVWell4-MseI      41     41    37
## 96       2           Roth-MLV3p-CD4TMLVwell4-MuA      21     21     1
## 97       2        Roth-MLV3p-CD4TMLVWell4-NlaIII      31     29    26
## 98       2       Roth-MLV3p-CD4TMLVWell4-Tsp509I      27     26    26
## 99       2         Roth-MLV3p-CD4TMLVWell5-BstYI      19     19    18
## 100      2          Roth-MLV3p-CD4TMLVWell5-MseI      40     40    37
## 101      2           Roth-MLV3p-CD4TMLVwell5-MuA      36     36    10
## 102      2        Roth-MLV3p-CD4TMLVWell5-NlaIII      23     21    19
## 103      2       Roth-MLV3p-CD4TMLVWell5-Tsp509I      54     53    48
## 104      2         Roth-MLV3p-CD4TMLVWell6-BstYI      14     14    12
## 105      2          Roth-MLV3p-CD4TMLVWell6-MseI      33     33    31
## 106      2           Roth-MLV3p-CD4TMLVwell6-MuA      15     14     1
## 107      2        Roth-MLV3p-CD4TMLVWell6-NlaIII      23     23    20
## 108      2       Roth-MLV3p-CD4TMLVWell6-Tsp509I      46     46    41
## 109      3   Ocwieja-HIV896-CD4TND365-InfectionI     422    422   400
## 110      3  Ocwieja-HIV896-CD4TND365-InfectionII     295    293   279
## 111      3 Ocwieja-HIV896-CD4TND365-InfectionIII     269    269   253
## 112      3      Ocwieja-HIV896-CD4TND365-NoVirus      NA     NA    NA
##     vectored linkered psl sites
## 1         NA       NA  NA    NA
## 2          1       15   5     2
## 3         NA       12  NA    NA
## 4         NA       NA  NA    NA
## 5         NA       NA  NA    NA
## 6         NA       NA  NA    NA
## 7         NA        3   1     1
## 8         NA        7   1     1
## 9         NA        9   2     2
## 10        NA        4   1     1
## 11        NA        6   1    NA
## 12        NA        2  NA    NA
## 13        NA        2  NA    NA
## 14        NA       NA  NA    NA
## 15        NA        3  NA    NA
## 16        NA        4  NA    NA
## 17        NA        3  NA    NA
## 18        NA        1  NA    NA
## 19        NA        1  NA    NA
## 20        NA       NA  NA    NA
## 21        NA       NA  NA    NA
## 22        NA       19   9     6
## 23        NA       26  31    24
## 24         4       36  38    29
## 25        NA       34  31    22
## 26        NA       32  24    16
## 27         1       42  49    34
## 28         4       31  37    25
## 29         1       39  41    32
## 30        NA       15  17    11
## 31        NA       18  23    18
## 32        NA       42  74    67
## 33        NA       35  46    39
## 34        NA       39  58    46
## 35        NA        9   7     6
## 36        NA        5   5     5
## 37        NA       13  14    10
## 38        NA        8  16    10
## 39        NA        1   1     1
## 40        NA        8   8     4
## 41         1       38  51    42
## 42        NA       17  NA    NA
## 43        NA       22  39    32
## 44        NA       30  44    34
## 45        NA        9  18    15
## 46        NA       43  56    48
## 47        NA       16   3     2
## 48         1       28  33    22
## 49        NA       67  80    60
## 50        NA       17  10     8
## 51        NA       36  44    33
## 52        NA       13  NA    NA
## 53        NA       22  31    25
## 54         1       54  63    44
## 55        NA       NA  NA    NA
## 56        NA       17   5     2
## 57        NA        2  NA    NA
## 58        NA       NA  NA    NA
## 59        NA       NA  NA    NA
## 60        NA       NA  NA    NA
## 61        NA        4   1     1
## 62        NA       15   1     1
## 63        NA        6   2     2
## 64        NA        4   1     1
## 65        NA        4   1    NA
## 66        NA       NA  NA    NA
## 67        NA       NA  NA    NA
## 68        NA       NA  NA    NA
## 69        NA        2  NA    NA
## 70        NA        7  NA    NA
## 71        NA        4  NA    NA
## 72        NA        1  NA    NA
## 73        NA       NA  NA    NA
## 74        NA       NA  NA    NA
## 75        NA       NA  NA    NA
## 76        NA       14   9     6
## 77         1       20  31    24
## 78        NA       28  38    29
## 79        NA       24  31    22
## 80        NA       45  24    16
## 81        NA       46  49    34
## 82         4       36  37    25
## 83        NA       33  41    32
## 84         1       16  17    11
## 85        NA       17  23    18
## 86        NA       64  74    67
## 87         1       32  46    39
## 88        NA       42  58    46
## 89        NA        3   7     6
## 90        NA        3   5     5
## 91        NA       13  14    10
## 92        NA        9  16    10
## 93        NA       NA   1     1
## 94         1       13   8     4
## 95        NA       40  51    42
## 96        NA       19  NA    NA
## 97        NA       30  39    32
## 98        NA       27  44    34
## 99        NA       19  18    15
## 100       NA       40  56    48
## 101       NA       27   3     2
## 102       NA       23  33    22
## 103        1       54  80    60
## 104       NA       14  10     8
## 105       NA       33  44    33
## 106       NA       15  NA    NA
## 107        1       23  31    25
## 108       NA       45  63    44
## 109       66      357 281   210
## 110       48      233 201   145
## 111       41      188 171   114
## 112       NA       NA  NA    NA
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
## R version 3.6.0 (2019-04-26)
## Platform: x86_64-apple-darwin18.5.0 (64-bit)
## Running under: macOS Mojave 10.14.4
## 
## Matrix products: default
## BLAS/LAPACK: /usr/local/Cellar/openblas/0.3.5/lib/libopenblasp-r0.3.5.dylib
## 
## locale:
## [1] C/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] hiReadsProcessor_1.19.1     hiAnnotator_1.18.0         
##  [3] GenomicAlignments_1.20.0    Rsamtools_2.0.0            
##  [5] SummarizedExperiment_1.14.0 DelayedArray_0.10.0        
##  [7] BiocParallel_1.17.19        matrixStats_0.54.0         
##  [9] Biobase_2.44.0              GenomicRanges_1.36.0       
## [11] GenomeInfoDb_1.20.0         Biostrings_2.52.0          
## [13] XVector_0.24.0              IRanges_2.18.0             
## [15] S4Vectors_0.22.0            BiocGenerics_0.30.0        
## 
## loaded via a namespace (and not attached):
##  [1] pkgload_1.0.2          splines_3.6.0          foreach_1.4.4         
##  [4] assertthat_0.2.1       BSgenome_1.52.0        GenomeInfoDbData_1.2.1
##  [7] cellranger_1.1.0       remotes_2.0.4          sessioninfo_1.1.1     
## [10] pillar_1.3.1           backports_1.1.4        lattice_0.20-38       
## [13] glue_1.3.1             digest_0.6.18          colorspace_1.4-1      
## [16] Matrix_1.2-17          plyr_1.8.4             XML_3.98-1.19         
## [19] pkgconfig_2.0.2        devtools_2.0.2         zlibbioc_1.30.0       
## [22] purrr_0.3.2            scales_1.0.0           processx_3.3.0        
## [25] tibble_2.1.1           ggplot2_3.1.1          usethis_1.5.0         
## [28] withr_2.1.2            lazyeval_0.2.2         cli_1.1.0             
## [31] magrittr_1.5           crayon_1.3.4           readxl_1.3.1          
## [34] evaluate_0.13          memoise_1.1.0          ps_1.3.0              
## [37] fs_1.3.0               xml2_1.2.0             pkgbuild_1.0.3        
## [40] tools_3.6.0            prettyunits_1.0.2      stringr_1.4.0         
## [43] munsell_0.5.0          callr_3.2.0            compiler_3.6.0        
## [46] rlang_0.3.4            grid_3.6.0             RCurl_1.95-4.12       
## [49] iterators_1.0.10       rstudioapi_0.10        bitops_1.0-6          
## [52] gtable_0.3.0           codetools_0.2-16       roxygen2_6.1.1        
## [55] R6_2.4.0               knitr_1.22             dplyr_0.8.0.1         
## [58] rtracklayer_1.44.0     sonicLength_1.4.6      commonmark_1.7        
## [61] rprojroot_1.3-2        desc_1.2.0             stringi_1.4.3         
## [64] Rcpp_1.0.1             xfun_0.6               tidyselect_0.2.5
```

