#########################################################################################
# R script to generate tuple scores/stats for verification
#
# BS-seq data set with 4 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 4 cancer lesions with normal tissue
# 
#
# Stephany Orjuela, July 2018
#########################################################################################

library(GenomicRanges)
library(SummarizedExperiment)

#load("data/derASM_fullCancer.RData")
load("data/derASM_fullcancer2.RData")
load("data/tupleASM_fullCancer.RData")

#### Get CpG and snp keys from derived-asm table ####--------------------------------------------------------
#load("trueASM.sampleTable.allChroms_fullcov.RData") old
verbose <- T

#Get keys for derived asm
keyGR <- rowRanges(derASM) #1,209,052 new: 1,225,491
start(keyGR) <- end(keyGR) <- start(keyGR) - 1


#Get keys for ASM score
asmscoreGR <- rowRanges(ASM) #3,753,430 new:3,589,472

#### FILTER: Get tuples from asm-score-table that overlap with sites from derived-score-table ####-----

asmscoreGR_end <- subsetByOverlaps(asmscoreGR, keyGR) #820,106 for these I'll draw the curve

#### Get the snp.key for the tuples ####-----------------------------------------------------------

w <- findOverlaps(keyGR, asmscoreGR_end)
wtable <- as(w, "data.frame")
wtable <- wtable[order(wtable$subjectHits),]
tuplesl <- unique(wtable$subjectHits)

snp <- assay(derASM, "snp.table")

snp.tuple.table <- mclapply(tuplesl, function(g){
  w <- as.integer(which(wtable$subjectHits == g))

  if(length(w) == 1){
    return(snp[wtable$queryHits[w],])
  }
  if(length(w) > 1){
    return(snp[wtable$queryHits[w],][1,])
  }
}, mc.cores = 6)
#)

x <- unlist(snp.tuple.table)
posed.snp.tuple.table <- matrix(x, nrow = length(tuplesl), ncol = 12, byrow = T)
colnames(posed.snp.tuple.table) <- colnames(derASM)
#posed.snp.tuple.table <- t(snp.tuple.table)

#### Get tuple coverage data ####-----------------------------------------------------------------

#IF you don't take this out, the loop takes AGES
ref.cov <- assay(derASM,"ref.cov")
alt.cov <- assay(derASM,"alt.cov")
ref.meth <- assay(derASM,"ref.meth")
alt.meth <- assay(derASM,"alt.meth")


tuple.ref.cov <- t(vapply(tuplesl, function(g){
  w <- as.integer(which(wtable$subjectHits == g))
  
  #Merge tuple coverage for read and alt reads
  if(length(w) == 1){
    return(ref.cov[wtable$queryHits[w],])
  } 
  if(length(w) > 1){
    return(colSums(ref.cov[wtable$queryHits[w],]))
  }
  
}, double(12)))

tuple.alt.cov <- t(vapply(tuplesl, function(g){
  w <- as.integer(which(wtable$subjectHits == g))
  
  #Merge tuple coverage for read and alt reads
  if(length(w) == 1){
    return(alt.cov[wtable$queryHits[w],])
  } 
  if(length(w) > 1){
    return(colSums(alt.cov[wtable$queryHits[w],]))
  }
  
}, double(12)))

tuple.ref.meth <- t(vapply(tuplesl, function(g){
  w <- as.integer(which(wtable$subjectHits == g))
  
  #Merge tuple coverage for read and alt reads
  if(length(w) == 1){
    return(ref.meth[wtable$queryHits[w],])
  } 
  if(length(w) > 1){
    return(colSums(ref.meth[wtable$queryHits[w],]))
  }
  
}, double(12)))

tuple.alt.meth <- t(vapply(tuplesl, function(g){
  w <- as.integer(which(wtable$subjectHits == g))
  
  #Merge tuple coverage for read and alt reads
  if(length(w) == 1){
    return(alt.meth[wtable$queryHits[w],])
  } 
  if(length(w) > 1){
    return(colSums(alt.meth[wtable$queryHits[w],]))
  }
  
}, double(12)))


#### Get the derived ASM for the same tuples ####-------------------------------------------
# #To input to get the derived-true asm

posed.true.tuple.asm <- abs((tuple.alt.meth / tuple.alt.cov) - (tuple.ref.meth / tuple.ref.cov))

#new version of trueASM
prop <- (tuple.ref.meth + tuple.alt.meth) / (tuple.ref.cov + tuple.alt.cov)

stat.tuple.asm  <- ((tuple.ref.meth/tuple.ref.cov) - (tuple.alt.meth/tuple.alt.cov)) /
  sqrt(prop * (1 - prop) * ((1/tuple.ref.cov) + (1/tuple.alt.cov)))

#### Build summ Experiment ####----------------------------------------------------------------

tuple.derASM <- SummarizedExperiment(assays=S4Vectors::SimpleList(der.ASM = posed.true.tuple.asm,
                                                                  stat.ASM = stat.tuple.asm,
                                                                  snp.table = posed.snp.tuple.table,
                                                                  ref.cov = tuple.ref.cov,
                                                                  alt.cov = tuple.alt.cov,
                                                                  ref.meth = tuple.ref.meth,
                                                                  alt.meth = tuple.alt.meth),
                                     rowRanges=asmscoreGR_end)

save(tuple.derASM, file = "data/tuplederASM_fullCancer2.RData")
