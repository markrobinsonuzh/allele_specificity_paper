#!/usr/bin/env Rscript

#########################################################################################
# R script to get DMRs using dmrseq (tBS-seq) 
#
# TBS-seq data CRCs Vs Norm
#
# Stephany Orjuela, August 2019
#########################################################################################

library(dmrseq)
library(BiSeq)
library(ggplot2)
library(DAMEfinder)

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

#### Compare DAMEs to DMRs ####
load("data/tupledames_cimp.RData")
load("data/tupledames_noncimp.RData")

#load("data/noncimpCRCsVsNORM_DMRs_dmrseq.RData")
#load("data/cimpCRCsVsNORM_DMRs_dmrseq.RData")

filt_over <- function(DMRsfilt, dames){
  #positive beta is hypo
  sum(DMRsfilt$beta < 0 & DMRsfilt$qval <= 0.05)
  sum(DMRsfilt$beta > 0 & DMRsfilt$qval <= 0.05) 
  DMRsfilthyper <- DMRsfilt[DMRsfilt$beta < 0 & DMRsfilt$qval <= 0.05]
  message("hyper DMRs:")
  print(length(DMRsfilthyper))
  DMRsfilthypo <- DMRsfilt[DMRsfilt$beta > 0 & DMRsfilt$qval <= 0.05]
  message("hypo DMRs:")
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

#Number for table 1 in paper
filt_over(DMRscimp, dames_cimp)
filt_over(DMRs, dames_noncimp)

#Plot mean methylation of CRC vs normal in top DAMEs over top DMRS
bsrel <- rawToRel(BSr) #26,959,049
dmrmeth <- methLevel(bsrel) 

bsrel.sub <- bsrel[rowSums(totalReads(BSr) >= 10) >= 10,] #2,856,397

grbs <- rowRanges(bsrel.sub)
dmrmeth <- methLevel(bsrel.sub)
group <- colData(BSr)$group

plot_methsPerreg <- function(crc, norm, regnum, file, DMRs, DAMEs){
  
  #Get average meth across groups per site
  dmrmethCIMP <- rowMeans(dmrmeth[,group == crc])
  dmrmethNORMCIMP <- rowMeans(dmrmeth[,group == norm])
  
  over <- findOverlaps(DMRs[1:regnum],grbs)
  cluster.ids <- 1:regnum
  
  #get average per DMR
  cimp_perreg <- sapply(cluster.ids, function(Index){ 
    mean(dmrmethCIMP[subjectHits(over)[queryHits(over) == Index]], na.rm = TRUE)
  })
  
  cimpnorm_perreg <- sapply(cluster.ids, function(Index){ 
    mean(dmrmethNORMCIMP[subjectHits(over)[queryHits(over) == Index]], na.rm = TRUE)
  })
  
  #get average per DAME
  
  gr_cimpdames <- GRanges(paste0("chr",DAMEs$chr),
                          IRanges(DAMEs$start, DAMEs$end))
  over <- findOverlaps(gr_cimpdames[1:regnum],grbs)
  
  cimp_perdame <- sapply(cluster.ids, function(Index){ 
    mean(dmrmethCIMP[subjectHits(over)[queryHits(over) == Index]], na.rm = TRUE)
  })
  
  cimpnorm_perdame <- sapply(cluster.ids, function(Index){ 
    mean(dmrmethNORMCIMP[subjectHits(over)[queryHits(over) == Index]], na.rm = TRUE)
  })
  
  
  # get average for overlapping regs (take DMR as the overlaped reg)
  overboth <- subsetByOverlaps(DMRs[1:regnum], gr_cimpdames[1:regnum])
  over <- findOverlaps(overboth,grbs)
  
  over_perdame <- sapply(cluster.ids, function(Index){ 
    mean(dmrmethCIMP[subjectHits(over)[queryHits(over) == Index]], na.rm = TRUE)
  })
  
  overnorm_perdame <- sapply(cluster.ids, function(Index){ 
    mean(dmrmethNORMCIMP[subjectHits(over)[queryHits(over) == Index]], na.rm = TRUE)
  })
  
  
  methtab <- data.frame(cimp = c(cimp_perreg,cimp_perdame, over_perdame),
                        norm = c(cimpnorm_perreg,cimpnorm_perdame, overnorm_perdame),
                        reg = c(rep("DMR",regnum),rep("DAME",regnum), 
                                rep("BOTH", length(overboth))))
  
  
  p <- ggplot(methtab) + 
    geom_point(aes(norm, cimp, color = reg), alpha = 0.5) + 
    coord_cartesian(xlim = 0:1, ylim = 0:1) +
    theme_bw()
  #ggsave(sprintf("curvesNscatters/%s",file))
  #return(p)
  return(methtab)
}

p1 <- plot_methsPerreg("cimp", "NORMAL.cimp", 1000, "methVals_DAMEvsDMR_cimp.png", 
                 DMRscimp, dames_cimp)

p2 <- plot_methsPerreg("non", "NORMAL.non", 1000, "methVals_DAMEvsDMR_non.png", 
                 DMRs, dames_noncimp)

m4 <- cowplot::plot_grid(p1,p2, ncol=1, nrow = 2, labels = c("CIMP","non-CIMP"))
ggsave("curvesNscatters/methVals_DAMEvsDMR_simes_both_1000.png", m4, width = 8, 
       height = 10)


#### get DAMEs that are not detected with DMRs ####
load("data/tupledames_cimp.RData")
sime_cimp <- GRanges(dames_cimp$chr, IRanges(dames_cimp$start, dames_cimp$end), 
                     pval = dames_cimp$pvalSimes,FDR = dames_cimp$FDR)
load("data/tupledames_noncimp.RData")
sime_non <- GRanges(dames_noncimp$chr, IRanges(dames_noncimp$start, dames_noncimp$end),
                    pval = dames_noncimp$pvalSimes,FDR = dames_noncimp$FDR)

load("data/tupledames_cimp_emp.RData")
load("data/tupledames_cimp_emp_repeat.RData")
emp_cimp <- GRanges(dames_cimp$chr, IRanges(dames_cimp$start, dames_cimp$end),
                    pval = dames_cimp$pvalEmp,FDR = dames_cimp$FDR)
load("data/tupledames_noncimp_emp.RData")
emp_non <- GRanges(dames_noncimp$chr, IRanges(dames_noncimp$start, dames_noncimp$end),
                   pval = dames_noncimp$pvalEmp,FDR = dames_noncimp$FDR)


DMRsfilt <- DMRs[DMRs$qval <= 0.05]
seqlevels(DMRsfilt) <- gsub("chr","",seqlevels(DMRsfilt))

sime_non <- sime_non[sime_non$FDR <= 0.05] #258

tes <- subsetByOverlaps(sime_non,DMRsfilt, invert = TRUE) #93

#Plot DAMEs not hitting DMRs
load("data/tupleASM_fullCancer.RData")
metadata <- read.table("data/fullCancerSamples.txt", stringsAsFactors = FALSE)
metadata$V2 <- ifelse(metadata$V2 == "NORM_non", "NORM_noncimp", metadata$V2)
metadata$V2 <-ifelse(metadata$V2 == "CRC_non", "CRC_noncimp", metadata$V2)
colData(ASM)$group <- metadata$V2  
colData(ASM)$samples <- colnames(ASM)
source("custom_scoretracks.R")

#non-CIMP
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(3,8)]
#9 136658255-136658387
a1 <- dame_track_forpap(dame = tes[2], ASM = ASM[,seq(1,11,2)], window = 2, colvec = myColor,
                        plotMeth = TRUE)
# b2 <- dame_track_forpap(dame = tes[3], ASM = ASM[,seq(1,11,2)], window = 3, colvec = myColor,
#                  plotMeth = TRUE)
# c3 <- dame_track(dame = tes[4], ASM = ASM[,seq(1,11,2)], window = 3, colvec = myColor)

emp_non <- emp_non[emp_non$FDR <= 0.05] #119 non

tes <- subsetByOverlaps(emp_non,DMRsfilt, invert = TRUE) #2

a2 <- dame_track_forpap(dame = tes[1], ASM = ASM[,seq(1,11,2)], window = 3, colvec = myColor,
                        plotMeth = TRUE)
# b2 <- dame_track_forpap(dame = tes[2], ASM = ASM[,seq(1,11,2)], window = 3, colvec = myColor,
#                  plotMeth =TRUE)

#### CIMP
DMRsfilt <- DMRscimp[DMRscimp$qval <= 0.05]
seqlevels(DMRsfilt) <- gsub("chr","",seqlevels(DMRsfilt))

sime_cimp <- sime_cimp[sime_cimp$FDR <= 0.05] #4051

tes <- subsetByOverlaps(sime_cimp,DMRsfilt, invert = TRUE)
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5)]
a3 <- dame_track_forpap(dame = tes[1], ASM = ASM[,seq(2,12,2)], window = 8, colvec = myColor,
                        plotMeth = TRUE)
b3 <- dame_track_forpap(dame = tes[3], ASM = ASM[,seq(2,12,2)], window = 1, colvec = myColor,
                        plotMeth = TRUE)
cowplot::plot_grid(a1,a2,a3,b3, nrow = 4, ncol = 1, align = "v", labels = c("A","","B",""))
ggplot2::ggsave(filename = "curvesNscatters/DAMEnotinDMRs_combi.png",
                width = 10, height = 12)


