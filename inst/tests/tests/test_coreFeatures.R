test_that("all core function tests", {

#### extractor tests ####
samplename <- "Roth-MLV3p-CD4TMLVWell6-Tsp509I"
subjectSeqs <- extractSeqs(seqProps, feature='genomicLinkered', sector = "1", 
                           samplename = samplename)[[1]][[1]]
subjectSeqs <- head(subjectSeqs)
expect_equal(as.character(subjectSeqs), 
             structure(c("CCTTCTACCACTGTCTTCAATT", "GGAGAATCCTTTTGTTATTTCCTGGAAATGCTTGAATCATAGATGAGTCTCTTCCTGCGTCCGTGTATTTCAGCTGCAGTTTCCAAGTGTGTAATAGTGCTGTCCACAAGGAATAGCTGCTGTTACTGTAAATT", "CTATGCCCTTCTGCGAGGTCAGAAAGGGCACGCACTCACACCCGTCCCACCTTCAATT", "ATATGAGTAGGACAGTTGGCAGATGAAATT", "TTATGAAGCCCGTTCTTCAGAAAGGTTTATAGCACTTTCGTTAAACCTCCTTCCTGTTCAAAATGCAGCGCACTAGAGTATTTTTCCCCAATCCACCGTTAAGCGGCTTCTTCTTAAGCAGGACTGTTTTGGGTCTCTTGGTTGCTACCCCTACTACACCCTACTAATATGTTGCACAAAGCGACGGGGAACCCCTTCTCCCCCAAACGATGTGCGTGGAAACACACAATGTCCCACTCTCTCCATTCCCAGAAACCTAGACTCTCCCCTCAAATACCTTTCCTGGCTCATCTACAGCACTAAGGTCTCAAGGTAGGGGCAAAAGCCCGCGACGATT", "GAACATCTCACCCGGATTTTTTACTGACTCCACTAAAATT"), .Names = c("GQKHUJB01AS3TK", "GQKHUJB01A92GD", "GQKHUJB01AXJRO", "GQKHUJB01A3C99", "GQKHUJB01BCP6X", "GQKHUJB01AL0CB")))

subjectSeqs <- extractFeature(seqProps, feature='decoded', sector = "1", 
                              samplename = samplename)[[1]][[1]]
subjectSeqs <- head(subjectSeqs)
expect_equal(as.character(subjectSeqs), 
             structure(c("TGATTGACTACCCGTCAGCGGGGGTCTTTCACCTTCTACCACTGTCTTCAATTGTCCCTTAAGCGGAGCCC", "TGATTGACTACCCGTCAGCGGGGGGTCTTTCAGGAGAATCCTTTTGTTATTTCCTGGAAATGCTTGAATCATAGATGAGTCTCTTCCTGCGTCCGTGTATTTCAGCTGCAGTTTCCAAGTGTGTAATAGTGCTGTCCACAAGGAATAGCTGCTGTTACTGTAAATTGTCCCTTAAGCGGAGCCC", "TGATTGACTACCCGTCAGCGGGGGTCTTTCACTATGCCCTTCTGCGAGGTCAGAAAGGGCACGCACTCACACCCGTCCCACCTTCAATTGTCCCTTAAGCGGAGCCC", "TGATTGACTACCCGTCAGCGGGGGTCTTTCAATATGAGTAGGACAGTTGGCAGATGAAATTGTCCCTTAAGCGGAGCCC", "TGATTGACTACCCGTCAGCGGGGGGTCTTTCATTATGAAGCCCGTTCTTCAGAAAGGTTTATAGCACTTTCGTTAAACCTCCTTCCTGTTCAAAATGCAGCGCACTAGAGTATTTTTCCCCAATCCACCGTTAAGCGGCTTCTTCTTAAGCAGGACTGTTTTGGGTCTCTTGGTTGCTACCCCTACTACACCCTACTAATATGTTGCACAAAGCGACGGGGAACCCCTTCTCCCCCAAACGATGTGCGTGGAAACACACAATGTCCCACTCTCTCCATTCCCAGAAACCTAGACTCTCCCCTCAAATACCTTTCCTGGCTCATCTACAGCACTAAGGTCTCAAGGTAGGGGCAAAAGCCCGCGACGATTGTCCCTTAAGCGGA", "TGATTGACTACCCGTCAGCGGGGGTCTTTCAGAACATCTCACCCGGATTTTTTACTGACTCCACTAAAATTGTCCCTTAAGCGGAGCCC"), .Names = c("GQKHUJB01AS3TK", "GQKHUJB01A92GD", "GQKHUJB01AXJRO", "GQKHUJB01A3C99", "GQKHUJB01BCP6X", "GQKHUJB01AL0CB")))

patternSeq <- extractFeature(seqProps, feature='primerltrsequence' ,sector = "1", 
                              samplename = samplename)
patternSeq <- patternSeq[[1]][[1]]
expect_equal(patternSeq, "TGATTGACTACCCGTCAGCGG")

#### test pairwise alignments to primers ####
res <- pairwiseAlignSeqs(subjectSeqs, patternSeq, qualityThreshold=0.85)
expect_that(res, is_identical_to(
  new("IRanges", start = c(1L, 1L, 1L, 1L, 1L, 1L) , 
      width = c(21L, 21L, 21L, 21L, 21L, 21L), 
      NAMES = c("GQKHUJB01AS3TK", "GQKHUJB01A92GD", "GQKHUJB01AXJRO", 
                "GQKHUJB01A3C99", "GQKHUJB01BCP6X", "GQKHUJB01AL0CB"), 
      elementType = "integer", elementMetadata = NULL, metadata = list())))

#### test pairwise alignments to linkers ####
patternSeq <- extractFeature(seqProps, feature='linkersequence' ,sector = "1", 
                             samplename = samplename)
patternSeq <- patternSeq[[1]][[1]]
expect_equal(patternSeq, "GTCCCTTAAGCGGAGCCCT")

res <- pairwiseAlignSeqs(subjectSeqs, patternSeq, qualityThreshold=0.55,
                         side="middle")
expect_that(res, is_identical_to(
  new("IRanges", 
      start = c(54L, 167L, 90L, 62L, 370L, 72L), 
      width = c(18L, 18L, 18L, 18L, 14L, 18L), 
      NAMES = c("GQKHUJB01AS3TK", "GQKHUJB01A92GD", "GQKHUJB01AXJRO", 
                "GQKHUJB01A3C99", "GQKHUJB01BCP6X", "GQKHUJB01AL0CB"),
      elementType = "integer", elementMetadata = NULL , metadata = list())))

#### test pairwise alignments to linkers with primerID ####
ids <- c("GGTTCTACGT", "AGGAGTATGA", "TGTCGGTATA", "GTTATAAAAC", "AGGCTATATC", 
         "ATGGTTTGTT")
subjectSeqs <- xscat(subjectSeqs, xscat(ids,"TTTTTTTTTTT"))
patternSeq <- "AAGCGGAGCCCNNNNNNNNNNTTTTTTTTTTT"
res <- primerIDAlignSeqs(subjectSeqs, patternSeq, doAnchored = TRUE)
expect_that(res, is_identical_to(
  new("CompressedIRangesList", elementType = "IRanges", elementMetadata = NULL,
      metadata = list(), 
      partitioning = new("PartitioningByEnd", end = c(5L, 10L), 
                         NAMES = c("hits", "primerIDs"), elementType = "integer"
                         , elementMetadata = NULL, metadata = list()), 
      unlistData = new("IRanges", 
                       start = c(61L, 174L, 97L, 69L, 79L, 72L, 185L, 108L, 
                                 80L, 90L), 
                       width = c(32L, 32L, 32L, 32L, 32L, 10L, 10L, 10L, 10L, 
                                 10L),
                       NAMES = c("read 1", "read 2", "read 3", "read 4", 
                                 "read 6", "read 1", "read 2", "read 3", 
                                 "read 4", "read 6"), elementType = "integer",
                       elementMetadata = NULL, metadata = list()
                      ))))

#### test writers & readers ####
psl <- head(extractFeature(seqProps, feature='psl', sector = "1", 
                           samplename = samplename)[[1]][[1]])
expect_that(write.psl(psl), prints_text("out.psl"))
expect_that(file.remove("out.psl"), is_true())

#### test utils ####
})