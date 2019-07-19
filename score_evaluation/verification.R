###########################################################################################
# R script to compare true snp-ASM to the tuple-ASM in a single sample. Make some roc 
# curves
#
# BS-seq data set with 6 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 6 cancer lesions with normal tissue
# 
#
# Stephany Orjuela, February 2018
###########################################################################################

library(GenomicRanges)
library(SummarizedExperiment)
library(iCOBRA)
library(ggplot2)

####Get ASM score ####---------------------------------------------------------------------

#load("tuple.trueASM.sampleTable.allChroms_fullcov.RData") #old
#load("trueASM.sampleTable.allChroms_fullcov.RData") #old
load("data/derASM_fullCancer.RData") #single site derived-true asm
load("data/tupleASM_fullCancer.RData") #tuple derived-true asm
load("data/tuplederASM_fullCancer.RData") #<-- from generateTupleTables.R script

keyGR <- rowRanges(derASM) #1,225,491
start(keyGR) <- end(keyGR) <- start(keyGR) - 1
rm(derASM)

asmscoreGR <- rowRanges(ASM) #3,589,472
over <- findOverlaps(asmscoreGR, keyGR)
ASM <- ASM[unique(sort(queryHits(over))),] #820,106

asmscoreGR <- rowRanges(ASM)

#read in amrfinder score
allelicfile <- "../amrfinder/NORM1.allelic"
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


#number = sample index
#above = vector of lower limit coverages
#below =  vector of upper limit coverages

x <- chooseSample(tuple.derASM, ASM, asmscoreGR, 7, 
                  c(5,10,50), c(9,49,3000), am)

#### Plot scores with eachother ####

xtab <- as.data.frame(x) #237,323
xtab$coverage <- gsub("50-3000", ">= 50", xtab$coverage)
xtab$coverage <- factor(xtab$coverage, 
                       levels = c("5-9","10-49", ">= 50"))
myColor <- RColorBrewer::brewer.pal(9, "Set1")

#asmtuple
p1 <- ggplot(xtab, aes(asm_tuple,derTrue_asm)) +
  geom_point(color = myColor[1]) +
  facet_wrap(~coverage) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines"))+
  labs(title = "ASMtuple", x = "ASMtuple", y = "ASMsnp")

#allelicmeth
p2 <- ggplot(xtab, aes(allelicmeth,derTrue_asm)) +
  geom_point(color = myColor[2]) +
  facet_wrap(~coverage) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines")) +
  labs(title = "allelicmeth", x = "allellicmeth", y = "ASMsnp")

#beta
p3 <-  ggplot(xtab, aes(beta,derTrue_asm)) +
  geom_point(color = myColor[3]) +
  facet_wrap(~coverage) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines")) +
  labs(title = "Scaled beta", x = "scaled beta", y = "ASMsnp")

p4 <- cowplot::plot_grid(p1,p2,p3, ncol=1, nrow = 3, align="v")
ggplot2::ggsave("curvesNscatters/scoreScatters.png", p4, width = 12, height = 12)

### Make curve ####-------------------------------------------------------------------------

#choose truth threshold (cases of asm)
#real1 <- ifelse(x$derTrue_asm >= 0.2, 1, 0) 
real2 <- ifelse(x$derTrue_asm >= 0.5, 1, 0)
real3 <- ifelse(x$derTrue_asm >= 0.8, 1, 0)

#generate truth + facet table
truth <- data.frame(#real1 = real1,
                    real2 = real2,
                    real3 = real3,
                    allreal = (real2+real3),
                    coverage = x$coverage,
                    ASM_snp = x$derTrue_asm)
                    #chromosome = seqnames(x))

#run iCOBRa for different true-thresholds
cobradat <- COBRAData(score = as.data.frame(mcols(x)[,3:5]),
                      truth = truth)

# cobraperf <- calculate_performance(cobradat, binary_truth = "real1", 
#                                    cont_truth = "ASM_snp", splv = "coverage",
#                                    aspects = "roc", maxsplit = Inf)
# 
# roc1 <- cobraperf@roc
# roc1$real <- "ASMsnp >= 0.2"

cobraperf <- calculate_performance(cobradat, binary_truth = "real2", 
                                   cont_truth = "ASM_snp", splv = "coverage",
                                   aspects = "roc", maxsplit = Inf)

roc2 <- cobraperf@roc
roc2$real <- "ASMsnp >= 0.5"

cobraperf <- calculate_performance(cobradat, binary_truth = "real3", 
                                   cont_truth = "ASM_snp", splv = "coverage",
                                   aspects = "roc", maxsplit = Inf)

roc3 <- cobraperf@roc
roc3$real <- "ASMsnp >= 0.8"

#put together all tables and plot
roc <- rbind(roc2,roc3)
roc <- roc[roc$splitval != "overall",]
roc$splitval <- gsub("coverage:5-9", "coverage = 5-9", roc$splitval)
roc$splitval <- gsub("coverage:10-49", "coverage = 10-49", roc$splitval)
roc$splitval <- gsub("coverage:50-3000", "coverage >= 50", roc$splitval)

roc$splitval <- factor(roc$splitval, 
                       levels = c("coverage = 5-9","coverage = 10-49", "coverage >= 50"))

roc$method <- gsub("asm_tuple", "ASMtuple",roc$method)
roc$method <- factor(roc$method, levels = c("ASMtuple","allelicmeth","beta"))

ggplot(roc, aes(FPR,TPR, color = method)) + 
  geom_line(size = 1) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines"), 
        text = element_text(size = 12)) +
  labs(color = "Score") +
  facet_wrap(~real+splitval, nrow = 3, ncol = 3) +
  scale_color_manual(values = myColor)
ggsave("curvesNscatters/full_ROCs_icobraggplot_reduced.png")
