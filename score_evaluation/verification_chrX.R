###########################################################################################
# R script to compare ASM_tuple in chrX Vs somatic genome. Make some ROC 
# curves
#
# WGBS-seq data set with 6 samples, 3 females and 3 males
# 
#
# Stephany Orjuela, September 2019
###########################################################################################

library(GenomicRanges)
library(SummarizedExperiment)
library(iCOBRA)
library(ggplot2)

####Get ASM score ####---------------------------------------------------------------------
load("data/malevsfem_ASM_mat.RData")
ASM_mat <- ASM_mat[rowSums(
  !is.na(assays(ASM_mat)[["cov"]]) &
    assays(ASM_mat)[["cov"]] >= 10) == BiocGenerics::ncol(ASM_mat),]

keyGR <- rowRanges(ASM_mat) #3,198,678

#read in allelicmeth score
allelicfile <- "../../../Shared_s3it/sorjuela/EGAdata/EGAD00001002523/C0010K/C0010K.allelic"
am <- data.table::fread(allelicfile, select = c(1,2,5))

#read in amrfinder score
amrfile <- "../../../Shared_s3it/sorjuela/EGAdata/EGAD00001002523/C0010K/C0010K.amr"
amr <- data.table::fread(amrfile, select = c(1:3,5))


#Choose a sample to compare the score
chooseSample <- function(ASM_score_matrix, scoreGR,
                         number, allelicfile, amrfile){
  
  
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
  
  
  print(length(scoreGR))
  
  #Remove weird Y chroms in females
  scoreGR <- scoreGR[seqnames(scoreGR) != "Y"]
  print(length(scoreGR))
  
  return(scoreGR)
}


#number = sample index

x <- chooseSample(ASM_mat, keyGR, 4,  am, amr)


### Make curve ####-------------------------------------------------------------------------

#choose truth threshold (in X chrom)
proms <- annotatr::build_annotations(genome = 'hg38', annotations = "hg38_genes_promoters")
seqlevels(proms) <- gsub("chr", "", seqlevels(proms))
xproms <- proms[seqnames(proms) == "X"]

over <- findOverlaps(x, xproms) #4,725
uniQ <- unique(queryHits(over))
real2 <- integer(length(x))
real2[uniQ] <- 1

#get AUCs
methodsl <- list(x$asm_tuple,x$beta, x$allelicmeth, x$amr)
auc1 <- unlist(lapply(methodsl, pROC::auc, response = real2))


#generate truth + facet table
truth <- data.frame(
  real2 = real2)

score <- as.data.frame(mcols(x)[,2:5])
rownames(truth) <- rownames(score)

#Plot precision(PPV)-recall(TPR)
library("ROCR")
truth2 <- data.frame(asm_tuple = truth$real2,
                     beta = truth$real2,
                     allelicmeth = truth$real2,
                     amr = truth$real2) 
pred <- prediction(score,truth2)
perf <- performance(pred,"prec","rec")

pr.dat <- data.frame(Recall = c(perf@x.values[[1]],perf@x.values[[2]],
                          perf@x.values[[3]],perf@x.values[[4]]),
           Score = factor(c(rep("ASMtuple",length(perf@x.values[[1]])), 
                     rep("methdeviation",length(perf@x.values[[2]])), 
                     rep("allelicmeth",length(perf@x.values[[3]])),
                     rep("amrfinder",length(perf@x.values[[4]]))),
                     levels = c("ASMtuple","methdeviation","allelicmeth","amrfinder")),
           Precision = c(perf@y.values[[1]],perf@y.values[[2]],
                             perf@y.values[[3]],perf@y.values[[4]]))


myColor <- RColorBrewer::brewer.pal(9, "Set1")
ggplot(pr.dat, aes(Recall,Precision, color = Score)) + 
  geom_line(size = 1) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines"), 
        text = element_text(size = 12)) +
  labs(color = "Score") +
  scale_color_manual(values = myColor)
ggsave("curvesNscatters/PrecRecall_chromX_proms.png")

#run iCOBRa for different true-thresholds
cobradat <- COBRAData(score = score,
                      truth = truth)

cobraperf <- calculate_performance(cobradat, binary_truth = "real2", 
                                   aspects = "roc")

roc <- cobraperf@roc


#put together all tables and plot

roc$method <- gsub("asm_tuple", "ASMtuple",roc$method)
roc$method <- gsub("amr", "amrfinder",roc$method)
roc$method <- gsub("beta", "methdeviation",roc$method)
roc$method <- factor(roc$method, levels = c("ASMtuple","methdeviation","allelicmeth","amrfinder"))


ggplot(roc, aes(FPR,TPR, color = method)) + 
  geom_line(size = 1) +
  theme_bw() + 
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        panel.spacing = unit(0, "lines"), 
        text = element_text(size = 12)) +
  labs(color = "Score") +
  scale_color_manual(values = myColor)
ggsave("curvesNscatters/ROCs_chromX_proms.png")
