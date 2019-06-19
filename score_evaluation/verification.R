#########################################################################################
# R script to compare true snp-ASM to the tuple-ASM in a single sample. Make some roc curves
#
# BS-seq data set with 6 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 6 cancer lesions with normal tissue
# 
#
# Stephany Orjuela, February 2018
#########################################################################################

library(GenomicRanges)
library(SummarizedExperiment)
library(iCOBRA)
library(pROC)
library(ggplot2)
library(cowplot)


####Get ASM score ####---------------------------------------------------------------------

#load("tuple.trueASM.sampleTable.allChroms_fullcov.RData") #old
#load("trueASM.sampleTable.allChroms_fullcov.RData") #old
load("../data/derASM_fullCancer.RData") #single site derived-true asm
load("../data/tupleASM_fullCancer.RData") #tuple derived-true asm
load("../data/tuplederASM_fullCancer.RData") #<-- from generateTupleTables.R script

keyGR <- rowRanges(derASM) #1,209,052 new: 1,225,491
start(keyGR) <- end(keyGR) <- start(keyGR) - 1
rm(derASM)

asmscoreGR <- rowRanges(ASM) #3,753,430 new:3,589,472
over <- findOverlaps(asmscoreGR, keyGR)
ASM <- ASM[unique(sort(queryHits(over))),] #820,106

asmscoreGR <- rowRanges(ASM)

#read in amrfinder score
allelicfile <- "amrfinder/NORM1.allelic"
am <- data.table::fread(allelicfile, select = c(1,2,5))

#Choose a sample to compare the score
chooseSample <- function(tuple.derived_ASM_matrix, ASM_score_matrix, scoreGR,
                         number, above, below, amrfile){
  
  
  #Put derASM on asmscore GRanges
  scoreGR$derTrue_asm <- 0
  scoreGR$derTrue_asm <- abs(assays(tuple.derived_ASM_matrix)[["der.ASM"]][,number])
  
  #add asmscore
  scoreGR$asm_tuple <- 0
  scoreGR$asm_tuple <- assays(ASM_score_matrix)[["asm"]][,number]
  scoreGR$asm_tuple <- ifelse(is.na(scoreGR$asm_tuple),0,scoreGR$asm_tuple)
  
  #add beta values as a score
  meth_level1 <- (assays(ASM_score_matrix)[["MM"]][,number] + 
                    assays(ASM_score_matrix)[["MU"]][,number] + 
                    assays(ASM_score_matrix)[["UM"]][,number]) / assays(ASM_score_matrix)[["cov"]][,number]
  
  scoreGR$beta <- ifelse(meth_level1 <= 0.5 , 
                          meth_level1/0.5, 
                          (1-meth_level1)/0.5)
  
  scoreGR$beta <- ifelse(is.na(scoreGR$beta), 0, scoreGR$beta)
  
  
  #add allelicmeth score from amrfinder
  #the score is a pvalue, allthough in the methpipe that is not explicitely mentioned
  amGR <- GRanges(amrfile$V1, IRanges(start=amrfile$V2, width = 1), almeth = amrfile$V5)
  o <- findOverlaps(amGR, scoreGR)
  scoreGR$allelicmeth <- 0
  #transform pvalue for comparison
  #scoreGR$almeth[subjectHits(o)] <- qnorm(1-(amGR$almeth[queryHits(o)]/2))
  scoreGR$allelicmeth[subjectHits(o)] <- -log10(amGR$almeth[queryHits(o)]) 
  
  #For cov
  scoreGR$coverage <- 0
  
  for(i in 1:length(above)){
    
    #Choose coverage
    covsum <- assay(tuple.derived_ASM_matrix,"ref.cov")[,number] + 
      assay(tuple.derived_ASM_matrix,"alt.cov")[,number]
    
    covfilt <- covsum >= above[i] &
      covsum <= below[i]  &
      !is.na(covsum)
    #add factor
    scoreGR$coverage[covfilt] <- sprintf("%i-%i", above[i], below[i])
    
  }
  print(length(scoreGR))
  
  #Remove weird cases of NA cov with SNP ID &
  #For the eval, remove what doesnt have a snp in this sample
  
  scoreGR <- scoreGR[!is.na(assay(tuple.derived_ASM_matrix,"der.ASM")[,number]) &
                     !is.na(assays(tuple.derived_ASM_matrix)[["snp.table"]][,number])]
  
  print(length(scoreGR))
  
  return(scoreGR)
}


x <- chooseSample(tuple.derASM, ASM, asmscoreGR, 7, 
                  c(2,5,10,50,100), c(4,9,49,99,3000), am)

### Make curve ####-----------------------------------------------------------------------------

#choose truth threshold (cases of asm)
real1 <- ifelse(x$derTrue_asm >= 0.4, 1, 0) 
real2 <- ifelse(x$derTrue_asm >= 0.6, 1, 0)
real3 <- ifelse(x$derTrue_asm >= 0.8, 1, 0)

#change upper name
cov <- ifelse(x$coverage == "100-3000", ">100", x$coverage)
cov <- factor(cov, levels = c("5-9","10-49","50-99",">100"))

#generate truth + facet table
truth <- data.frame(real1 = real1,
                    real2 = real2,
                    real3 = real3,
                    allreal = (real1+real2+real3),
                    coverage = cov,
                    ASM_snp = x$derTrue_asm)
                    #chromosome = seqnames(x))

#run iCOBRa
cobradat <- COBRAData(score = as.data.frame(mcols(x)[,3:5]),
                      truth = truth)

#real2, by cov
cobraperf <- calculate_performance(cobradat, binary_truth = "real2", 
                                   cont_truth = "ASM_snp", splv = "coverage",
                                   aspects = "roc", maxsplit = Inf)

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Set1", 
                                   facetted = TRUE)

plot_roc(cobraplot, title = "ASMsnp = 0.6") + facet_wrap(~splitval, nrow = 1)

