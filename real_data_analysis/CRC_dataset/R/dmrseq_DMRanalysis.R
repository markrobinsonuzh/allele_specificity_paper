#!/usr/bin/env Rscript

#########################################################################################
# R script to get DMRs using dmrseq (tBS-seq) 
#
# TBS-seq data CRCs Vs Norm
#
# Stephany Orjuela, August 2019
#########################################################################################

library(dmrseq)

load("/home/sorjuela/serrated_pathway_paper/BiSulf/Data/CRCdata/BSrawCRC.RData")
bs <- BSseq(chr = seqnames(BSr), 
            pos = start(BSr),
            M = methReads(BSr), 
            Cov = totalReads(BSr), 
            sampleNames = colnames(BSr))

pData(bs)$group <- colData(BSr)$group
pData(bs)$patient <- rep(1:6, each = 2)
bs <- sort(bs)

#for some reason the old bsseq obj is more heavy than it should be
#so I re-create the BSseq object 

# library(bsseq)
# bs_new <- BSseq(M = getCoverage(bs, type = "M"),
#                 Cov = getCoverage(bs, type = "Cov"),
#                 pos = start(bs),
#                 chr = as.character(seqnames(bs)))
# pData(bs_new) <- pData(bs)

#I don't know how to specify contrasts in this thing, so i have to split the bs 
#object in two and run for each one

#CIMP comparison
bscimp <- bs[,pData(bs)$group %in% c("cimp","NORMAL.cimp")]
pData(bscimp)$group <- factor(pData(bscimp)$group, levels = c("cimp", "NORMAL.cimp"))

loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bscimp, type="Cov")==0) == 0)
bscimp <- bscimp[loci.idx,]

DMRscimp <- dmrseq(bs=bscimp, testCovariate="group", cutoff = 0.05, 
               BPPARAM = MulticoreParam(1),
               adjustCovariate = "patient")

save(DMRscimp, file = "data/cimpCRCsVsNORM_DMRs_dmrseq.RData")

#nonCIMP comparison
bsnon <- bs[,pData(bs)$group %in% c("non","NORMAL.non")]
pData(bsnon)$group <- factor(pData(bsnon)$group, levels = c("non", "NORMAL.non"))

loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bsnon, type="Cov")==0) == 0)
bsnon <- bsnon[loci.idx,]

DMRsnon <- dmrseq(bs=bsnon, testCovariate="group", cutoff = 0.05, 
                   BPPARAM = MulticoreParam(1),
                   adjustCovariate = "patient")

save(DMRsnon, file = "data/noncimpCRCsVsNORM_DMRs_dmrseq.RData")

#Load dames
load("data/tupledames_cimp.RData")
load("data/tupledames_noncimp.RData")

#Filter
#positive beta is hypo
filt_over <- function(DMRsfilt, dames){
  sum(DMRsfilt$beta < 0 & DMRsfilt$qval < 0.05)
  sum(DMRsfilt$beta > 0 & DMRsfilt$qval < 0.05) 
  DMRsfilthyper <- DMRsfilt[DMRsfilt$beta < 0 & DMRsfilt$qval < 0.05]
  print(length(DMRsfilthyper))
  DMRsfilthypo <- DMRsfilt[DMRsfilt$beta > 0 & DMRsfilt$qval < 0.05]
  print(length(DMRsfilthypo))
  
  GRcimp <- GRanges(dames$chr, 
                       IRanges(dames$start, dames$end),
                       A = dames$meanTstat,
                       pval = dames$pvalSimes,
                       FDR = dames$FDR)
  
  GRcimp <- GRcimp[GRcimp$FDR <= 0.05] #4037 cimp
  print(length(GRcimp))
  seqlevels(GRcimp) <- paste0("chr",seqlevels(GRcimp))
  over <- findOverlaps(GRcimp, DMRsfilthyper)
  over2 <- findOverlaps(GRcimp, DMRsfilthypo)
  
  message("dames that hit hyper DMRs ",length(GRcimp[unique(queryHits(over))]))
  
  message("Gain of ASM ",sum(GRcimp[unique(queryHits(over))]$A > 0))
  message("Loss of ASM ",sum(GRcimp[unique(queryHits(over))]$A < 0))
  
  message("dames that hit hypo DMRs ",length(GRcimp[unique(queryHits(over2))]))
  
  message("Gain of ASM ",sum(GRcimp[unique(queryHits(over2))]$A > 0))
  message("Loss of ASM ",sum(GRcimp[unique(queryHits(over2))]$A < 0))
  
  message("hyper DMRs that hit dames ",length(DMRsfilthyper[unique(subjectHits(over))]))
  message("hypo DMRs that hit dames ",length(DMRsfilthypo[unique(subjectHits(over2))]))
  
}
 
filt_over(DMRscimp, dames_cimp)
filt_over(DMRs, dames_noncimp)
