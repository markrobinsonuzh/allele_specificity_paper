#!/usr/bin/env Rscript

#########################################################################################
# R script plot several scores per region
# TBS-seq data CRCs Vs Norm
#
# Stephany Orjuela, September 2019
#########################################################################################

library(DAMEfinder)
library(ggplot2)
library(SummarizedExperiment)
library(cowplot)

metadata <- read.table("data/fullCancerSamples.txt", stringsAsFactors = FALSE)
metadata$V2 <- ifelse(metadata$V2 == "NORM_non", "NORM_noncimp", metadata$V2)
metadata$V2 <-ifelse(metadata$V2 == "CRC_non", "CRC_noncimp", metadata$V2)

#snp asm
load("data/derASM_fullcancer2.RData")
colData(derASM)$group <- metadata$V2  
colnames(derASM) <- c(paste0("C",1:6),paste0("N",1:6))
colData(derASM)$samples <- colnames(derASM)
colData(derASM)$samples == colnames(derASM)

#tuple asm
load("data/tupleASM_fullCancer.RData")
colData(ASM)$group <- metadata$V2  
colnames(ASM) <- c(paste0("C",1:6),paste0("N",1:6))
colData(ASM)$samples <- colnames(ASM)
colData(ASM)$samples == colnames(ASM)

#supp fig of opposite DAMEs
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5,1,3,4,8:9)]

DAME <- GRanges(9, IRanges(99984206,99984364))

a <- dame_track(dame = DAME, 
                window = 10, 
                derASM = derASM[,seq(2,12,2)],
                ASM = ASM[,seq(2,12,2)],
                colvec = myColor,
                plotSNP = TRUE
           ) 

ggplot2::ggsave(filename = "curvesNscatters/ASMopposDAME_func.png", plot = a,
                width = 8, height = 8)

snp <- GRanges(9, IRanges(99984349, width = 1))
allps <- mapply(methyl_circle_plot,
                vcfFile = gsub("/home/",path, metadata$V4)[seq(2,12,2)], 
                bamFile = gsub("/home/",path,metadata$V3)[seq(2,12,2)],
                sampleName = c(paste0("C",c(2,4,6)), paste0("N",c(2,4,6))),
                MoreArgs=list(
                  snp = snp,
                  refFile = gsub("/home/",path,reference_file),
                  dame = DAME,
                  sampleReads = TRUE,
                  numReads = 15,
                  pointSize = 2,
                  letterSize = 2.5
                ))

b <- cowplot::plot_grid(plotlist = allps, nrow = 2, ncol = 3)
ggplot2::ggsave(filename = "curvesNscatters/ASMopposDAME_reads.png", plot = b,
                width = 10, height = 8)

#figure 6
filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) >= 10
derASM <- derASM[filt,]

filt <- c(rowSums(assay(ASM, "cov") >= 10 & !is.na(assay(ASM, "cov"))) >= 10)
ASM <- ASM[filt,]

#MDS
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,3,5,8)]

m1 <- methyl_MDS_plot(derASM, group = metadata$V2, pointSize = 6, adj = 0.04) +
  scale_color_manual(values = myColor) +
  labs(color = "Tissue")

m2 <- methyl_MDS_plot(ASM, group = metadata$V2, pointSize = 6, adj = 0.05) +
  theme(legend.position = "none") + 
  scale_color_manual(values = myColor) 


#track
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5,1,3,4,8:9)]
DAME <- GRanges(9, IRanges(99983697,99984022))

## Add pvalue to the track, have to run limma first
grp <- factor(metadata$V2)
grp <- relevel(grp, "NORM_cimp")
samp <- gsub("CRC|NORM","", metadata$V1)
mod <- model.matrix(~0+grp)
mod <- mod[,-(7:9)] 


cont <- limma::makeContrasts(grpCRC_cimp-grpNORM_cimp, grpCRC_noncimp-grpNORM_noncimp, 
                             levels = mod)

# x <- assay(derASM, "der.ASM")
# fit <- limma::lmFit(x,mod)
# fit$coefficients <- fit$coefficients[,-(5:8)]
# fit$stdev.unscaled <- fit$stdev.unscaled[,-(5:8)]
# fit$cov.coefficients <- fit$cov.coefficients[,-(5:8)]
# fit$design <- fit$design[,-(5:8)] 
# 
# fit.cont <- limma::contrasts.fit(fit, contrasts = cont[-(5:8),])
# fit2 <- limma::eBayes(fit.cont)


tstats <- get_tstats(derASM, mod, contrast = cont, coef = 1, maxGap = 100, filter = FALSE)
tstatstup <- get_tstats(ASM, mod, contrast = cont, coef = 1, maxGap = 200, filter = FALSE)
source("custom_scoretracks.R")

m3pval <- dame_track_forpap(DAME, 
                  positions = 200,
                  derASM = tstats,
                  ASM = tstatstup
                  ) 
#m3pval
#ggplot2::ggsave("curvesNscatters/pval_track_fig6.png", m3pval, width=10, height = 4)


m3 <- dame_track(DAME, 
           #window = 5, 
           positions = 200,
           derASM = derASM[,seq(2,12,2)],
           ASM = ASM[,seq(2,12,2)],
           colvec = myColor,
           plotSNP = FALSE) +
  labs(color = "Tissue")

ggdraw() +
  draw_plot(m2, x = 0, y = .7, width = .42, height = .3) +
  draw_plot(m1, x = .42, y = .7, width = .58, height = .3) +
  draw_plot(m3, x = 0, y = 0.2,width = 1, height = 0.5) +
  draw_plot(m3pval, x = 0.021, y = 0,width = 0.81, height = 0.2) +
  draw_plot_label(label = c("A", "B", "C","D"), size = 13,
                  x = c(0, 0.42, 0,0), y = c(1, 1, 0.7,0.2))

ggplot2::ggsave("curvesNscatters/MDSboth_newfunc_withpval.png", width = 8, height = 9)

#figure 6
source("custom_mediantracks.R")

#MEG3
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5)]
DAME <- GRanges(14, IRanges(101291540,101293480))
megcimp <- dame_track_median_forpap(DAME, 
                 window = 10, 
                 ASM = ASM[,seq(2,12,2)],
                 colvec = myColor)

myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(3,8)]
megnon <- dame_track_median_forpap(DAME, 
                  window = 10, 
                  ASM = ASM[,seq(1,11,2)],
                  colvec = myColor)

#h19
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5)]
DAME <- GRanges(11, IRanges(2021017,2021260))
h19cimp <- dame_track_median_forpap(DAME, 
                      window = 10,
                      ASM = ASM[,seq(2,12,2)],
                      colvec = myColor)

myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(3,8)]
h19non <- dame_track_median_forpap(DAME, 
                                    window = 10, 
                                    ASM = ASM[,seq(1,11,2)],
                                    colvec = myColor)

#gnas
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5)]
DAME <- GRanges(20, IRanges(57425758,57428036))
gnascimp <- dame_track_median_forpap(DAME, 
                      positions = 700,
                      ASM = ASM[,seq(2,12,2)],
                      colvec = myColor)

myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(3,8)]
gnasnon <- dame_track_median_forpap(DAME, 
                       positions = 700,
                       ASM = ASM[,seq(1,11,2)],
                       colvec = myColor)

cowplot::plot_grid(megcimp, megnon, h19cimp, h19non, gnascimp, gnasnon, nrow = 6, 
                   ncol = 1, labels = c("A","","B","","C",""), align = "v")
ggplot2::ggsave(filename = "curvesNscatters/LOI_func_medians.png",
                width = 10, height = 12)

