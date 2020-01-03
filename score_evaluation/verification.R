###########################################################################################
# R script to compare true snp-ASM to the tuple-ASM per sample. Make some ROC 
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
load("data/derASM_fullcancer2.RData") #single site derived-true asm
load("data/tupleASM_fullCancer.RData") #tuple derived-true asm
load("data/tuplederASM_fullCancer2.RData") #<-- from generateTupleTables.R script

keyGR <- rowRanges(derASM) #1,453,873
start(keyGR) <- end(keyGR) <- start(keyGR) - 1
rm(derASM)

asmscoreGR <- rowRanges(ASM) #3,589,472
over <- findOverlaps(asmscoreGR, keyGR)
ASM <- ASM[unique(sort(queryHits(over))),] #820,106

asmscoreGR <- rowRanges(ASM)

#Already knowing which sample I'm going to plot
#read in allelicmeth score
allelicfile <- "../amrfinder/NORM1.allelic"
am <- data.table::fread(allelicfile, select = c(1,2,5))

#read in amrfinder score
amrfile <- "../amrfinder/NORM1.amr"
amr <- data.table::fread(amrfile, select = c(1:3,5))


#Choose a sample to compare the score
chooseSample <- function(tuple.derived_ASM_matrix, ASM_score_matrix, scoreGR,
                         number, above, below, allelicfile, amrfile){
  
  
  #Put derASM on asmscore GRanges
  scoreGR$derTrue_asm <- 0
  scoreGR$derTrue_asm <- abs(assays(tuple.derived_ASM_matrix)[["der.ASM"]][,number])
  
  scoreGR$stat_asm <- 0
  scoreGR$stat_asm <- abs(assays(tuple.derived_ASM_matrix)[["stat.ASM"]][,number])
  
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
  amGR <- GRanges(allelicfile$V1, IRanges(start=allelicfile$V2, width = 1), almeth = allelicfile$V5)
  start(amGR) <- start(amGR) + 1
  end(amGR) <- end(amGR) + 2
  o <- findOverlaps(amGR, scoreGR)
  scoreGR$allelicmeth <- 0
  scoreGR$allelicmeth[subjectHits(o)] <- -log10(amGR$almeth[queryHits(o)]) 
  scoreGR$allelicmeth <- ifelse(is.infinite(scoreGR$allelicmeth), 313.8136, scoreGR$allelicmeth)
  
  #add amrfinder score
  amrGR <- GRanges(amrfile$V1, IRanges(start=amrfile$V2, end=amrfile$V3), almeth = amrfile$V5)
  start(amrGR) <- start(amrGR) + 1
  end(amrGR) <- end(amrGR) + 2
  o <- findOverlaps(amrGR, scoreGR)
  scoreGR$amr <- 0
  scoreGR$amr[subjectHits(o)] <- -log10(amrGR$almeth[queryHits(o)]) 
  scoreGR$amr <- ifelse(is.infinite(scoreGR$amr), 313.8136, scoreGR$amr)
  
  
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
  
  #remove cov below 5
  scoreGR <- scoreGR[scoreGR$coverage != "0"]
  print(length(scoreGR))
  
  return(scoreGR)
}


#number = sample index
#above = vector of lower limit coverages
#below =  vector of upper limit coverages

x <- chooseSample(tuple.derASM, ASM, asmscoreGR, 7, #NORM1_non
                  c(5,10,50), c(9,49,4000), am, amr)

#### Plot scores distributions ####

xtab <- as.data.frame(x) #237,323, 282,097
xtab$coverage <- gsub("50-4000", ">= 50", xtab$coverage)
xtab$coverage <- factor(xtab$coverage, 
                        levels = c("5-9","10-49", ">= 50"))
myColor <- RColorBrewer::brewer.pal(9, "Set1")

# over <- findOverlaps(x,slopGR) #from ImprintedGenes script
# xtab_impr <- xtab[unique(queryHits(over)),]

## Distributions
xtabsub <- xtab[,c(7,9:13)]
colnames(xtabsub) <- c("ASMsnp","ASMtuple","methdeviation","allelicmeth","amrfinder",
                       "coverage")

xtabsub$threshold <- ifelse(xtabsub$ASMsnp >= 0.5, "ASMsnp >= 0.5", "ASMsnp < 0.5")
xtabsub <- reshape2::melt(xtabsub, 
                          id.vars = c("coverage", "threshold"),
                          variable.name = "Score")


p1 <- ggplot(xtabsub, aes(value, color = coverage, fill = coverage)) + 
  theme_bw() + 
  geom_density(alpha = 0.2) +
  theme(strip.background = element_rect(colour = "black", fill = "white")) +
        #panel.spacing = unit(0.5, "lines"), 
        #text = element_text(size = 12)) +
  facet_wrap(Score~threshold, scales = "free", nrow = 5) +
  scale_color_manual(values = c(myColor[5],myColor[1:4])) +
  scale_fill_manual(values = c(myColor[5],myColor[1:4])) +
  scale_x_continuous(trans='sqrt')
# threshold
# coverage ASMsnp < 0.5 ASMsnp >= 0.5
# 5-9           5909            67
# 10-49       108618           591
# >= 50       166286           626

xtabsub <- xtab[,c(7,9:13)]
colnames(xtabsub) <- c("ASMsnp","ASMtuple","methdeviation","allelicmeth","amrfinder",
                       "coverage")

#change threshold
xtabsub$threshold <- ifelse(xtabsub$ASMsnp >= 0.8, "ASMsnp >= 0.8", "ASMsnp < 0.8")
xtabsub <- reshape2::melt(xtabsub, 
                          id.vars = c("coverage", "threshold"),
                          variable.name = "Score")

p2 <- ggplot(xtabsub, aes(value, color = coverage, fill = coverage)) + 
  theme_bw() + 
  geom_density(alpha = 0.2) +
  theme(strip.background = element_rect(colour = "black", fill = "white")) +
  facet_wrap(Score~threshold, scales = "free", nrow = 5) +
  scale_color_manual(values = c(myColor[5],myColor[1:4])) +
  scale_fill_manual(values = c(myColor[5],myColor[1:4])) +
  scale_x_continuous(trans='sqrt')
# threshold
# coverage ASMsnp < 0.8 ASMsnp >= 0.8
# 5-9           5972             4
# 10-49       109169            40
# >= 50       166887            25

cowplot::plot_grid(p1,p2, labels = c("A","B"))
ggsave("curvesNscatters/scores_distributions_facet2.png",
       width = 10, height = 10)



## Scatters
#asmtuple
p1 <- ggplot(xtab, aes(asm_tuple,derTrue_asm)) +
  geom_point(color = myColor[1]) +
  geom_point(data=xtab_impr, aes(asm_tuple,derTrue_asm), color = myColor[6]) +
  scale_x_continuous(trans='sqrt') +
  scale_y_continuous(trans='sqrt') +
  facet_wrap(~coverage) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines"))+
  labs(title = "ASMtuple", x = "ASMtuple", y = "ASMsnp")

#allelicmeth
p2 <- ggplot(xtab, aes(allelicmeth,derTrue_asm)) +
  geom_point(color = myColor[3]) +
  geom_point(data=xtab_impr, aes(allelicmeth,derTrue_asm), color = myColor[6]) +
  scale_x_continuous(trans='sqrt') +
  scale_y_continuous(trans='sqrt') +
  facet_wrap(~coverage) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines")) +
  labs(title = "allelicmeth", x = "allellicmeth", y = "ASMsnp")

#beta
p3 <-  ggplot(xtab, aes(beta,derTrue_asm)) +
  geom_point(color = myColor[2]) +
  geom_point(data=xtab_impr, aes(beta,derTrue_asm), color = myColor[6]) +
  scale_x_continuous(trans='sqrt') +
  scale_y_continuous(trans='sqrt') +
  facet_wrap(~coverage) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines")) +
  labs(title = "Scaled beta", x = "scaled beta", y = "ASMsnp")

p4 <- cowplot::plot_grid(p1,p2,p3, ncol=1, nrow = 3, align="v")
ggplot2::ggsave("curvesNscatters/scoreScatters_withimpr.png", p4, width = 12, height = 12)

### Make curve ####-------------------------------------------------------------------------

#choose truth threshold (cases of asm)
real2 <- ifelse(x$derTrue_asm >= 0.5, 1, 0)
real3 <- ifelse(x$derTrue_asm >= 0.8, 1, 0)

#use new statASM
#real2 <- ifelse(x$stat_asm >= 7, 1, 0)
#real3 <- ifelse(x$stat_asm >= 10, 1, 0)

#get AUCs
methodsl <- list(x$asm_tuple,x$beta, x$allelicmeth, x$amr)
auc1 <- unlist(lapply(methodsl, pROC::auc, response = real2))
auc2 <- unlist(lapply(methodsl, pROC::auc, response = real3))

#generate truth + facet table
truth <- data.frame(
  real2 = real2,
  real3 = real3,
  allreal = (real2+real3),
  coverage = x$coverage,
  ASM_snp = x$derTrue_asm)

#run iCOBRa for different true-thresholds
cobradat <- COBRAData(score = as.data.frame(mcols(x)[,4:7]),
                      truth = truth)

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
roc$splitval <- gsub("coverage:50-4000", "coverage >= 50", roc$splitval)

roc$splitval <- factor(roc$splitval, 
                       levels = c("coverage = 5-9","coverage = 10-49", "coverage >= 50"))

roc$method <- gsub("asm_tuple", "ASMtuple",roc$method)
roc$method <- gsub("amr", "amrfinder",roc$method)
roc$method <- gsub("beta", "methdeviation",roc$method)
roc$method <- factor(roc$method, levels = c("ASMtuple","methdeviation","allelicmeth","amrfinder"))

ggplot(roc, aes(FPR,TPR, color = method)) + 
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = seq(0,0.75, 0.25)) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines"), 
        text = element_text(size = 12)) +
  labs(color = "Score") +
  facet_wrap(~real+splitval, nrow = 3, ncol = 3) +
  scale_color_manual(values = myColor)
ggsave("curvesNscatters/full_ROCs_icobraggplot_reduced_withamr_fix_methdev.png",
       width = 8, height = 4.5)
