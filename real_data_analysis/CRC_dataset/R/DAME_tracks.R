#September 2019

#How to explore DAMEs using DAMEfinder :)

library(GenomicRanges)
library(SummarizedExperiment)
library(cowplot)
#library(DAMEfinder)
devtools::load_all("../DAMEfinder_git/")

#Make sure this are done
load("data/tupledames_cimp.RData")
load("data/tupledames_noncimp.RData")
load("data/tupleASM_fullCancer.RData")
load("data/derASM_fullcancer2.RData")
metadata <- read.table("data/fullCancerSamples.txt", stringsAsFactors = FALSE)
colData(ASM)$group <- metadata$V2  
colData(derASM)$group <- metadata$V2  
colData(ASM)$samples <- colnames(ASM)

#I have a list of DAMEs like so, and I know what score I used to get it (e.g. ASMsnp or ASMtuple):
head(dames_cimp) #<- ASM_tuple

#inputs
DAME <- GRanges(9, IRanges(99984206,99984364))
ASM
derASM
myColor <- RColorBrewer::brewer.pal(9, "Set1")

#Step1
#Case 1.A plot with ASMtuple object
dame_track(DAME, ASM = ASM) + scale_color_manual(values = myColor)

#Case 1.B (if I have it) plot with ASMsnp object
dame_track(DAME, window = 10, derASM = derASM) + scale_color_manual(values = myColor)

#Case 1.C (If I have them) plot with both objects
dame_track(DAME, window = 10, derASM = derASM, ASM = ASM) + scale_color_manual(values = myColor)

#Step 2
reference_file <- "/home/Shared_taupo/data/annotation/Human/Ensembl_GRCh37.91/GRCh37.91.fa"
metadata <- read.table("data/fullCancerSamples.txt", stringsAsFactors = FALSE)
bam_files <- metadata$V3
vcf_files <- metadata$V4 
sample_names <- metadata$V1
path <- "/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/"

#Case 2.A
methyl_circle_plotCpG(cpgsite = DAME, #Here I would put the location of a CpG if i wanted
                      bamFile = gsub("/home/",path,metadata$V3[2]), 
                      refFile = gsub("/home/",path,reference_file),
                      dame = DAME, 
                      pointSize = 1,
                      order = TRUE)


methyl_circle_plotCpG(cpgsite = DAME, #Here I would put the location of a CpG if i wanted
                      bamFile = gsub("/home/",path,metadata$V3[2]), 
                      refFile = gsub("/home/",path,reference_file),
                      dame = DAME, 
                      pointSize = 1.5,
                      sampleName = "C2",
                      sampleReads = TRUE,
                      order = TRUE)


#Case 2.B 
snp <- GRanges(9, IRanges(99984349,width = 1))
methyl_circle_plot(snp, 
                   vcfFile = gsub("/home/",path, metadata$V4[2]),
                   bamFile = gsub("/home/",path,metadata$V3[2]), 
                   refFile = gsub("/home/",path,reference_file),
                   sampleName = "C2",
                   dame = DAME,
                   pointSize = 2,
                   sampleReads = TRUE)

methyl_circle_plot(snp, 
                   vcfFile = gsub("/home/",path, metadata$V4[8]),
                   bamFile = gsub("/home/",path,metadata$V3[8]), 
                   refFile = gsub("/home/",path,reference_file),
                   sampleName = "N2",
                   dame = DAME,
                   pointSize = 2,
                   sampleReads = TRUE)

##Do you think I should make a function to include all these plots at once?
## or provide code to do it? it a loop trough all the samples basically

#Case 2.C

allps <- mapply(methyl_circle_plot,
                vcfFile = gsub("/home/",path, metadata$V4), 
                bamFile = gsub("/home/",path,metadata$V3),
                sampleName = sample_names,
                MoreArgs=list(
                  snp = snp,
                  refFile = gsub("/home/",path,reference_file),
                  dame = DAME,
                  sampleReads = TRUE,
                  numReads = 15,
                  pointSize = 1,
                  letterSize = 1
                ))


cowplot::plot_grid(plotlist = allps, nrow = 2, ncol = 6)

# pdf("All_Mplots.pdf")
# lapply(allps, function(i) i)
# dev.off()


#### For paper ####

#supp fig of opposite DAMEs
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5,1,3,4,8:9)]

a <- dame_track(DAME, window = 10, 
           derASM = derASM[,seq(2,12,2)],
           ASM = ASM[,seq(2,12,2)],
           colvec = myColor) 

ggplot2::ggsave(filename = "curvesNscatters/ASMopposDAME_func.png", plot = a,
                width = 8, height = 8)

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

#figure 5
#MDS
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,3,5,8)]
colnames(derASM) <- c(paste0("C",1:6),paste0("N",1:6))
colnames(ASM) <- c(paste0("C",1:6),paste0("N",1:6))
colData(ASM)$samples <- colnames(ASM)
colData(derASM)$samples == colnames(derASM)

m1 <- methyl_MDS_plot(derASM, group = metadata$V2, pointSize = 6, adj = 0.04) +
  #theme(legend.position = "") + 
  scale_color_manual(values = myColor) +
  labs(color = "Tissue")

m2 <- methyl_MDS_plot(ASM, group = metadata$V2, pointSize = 6, adj = 0.05) +
  theme(legend.position = "none") + 
  scale_color_manual(values = myColor) 


#track
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5,1,3,4,8:9)]
DAME <- GRanges(9, IRanges(99983697,99984022))
m3 <- dame_track(DAME, 
           #window = 5, 
           positions = 200,
           derASM = derASM[,seq(2,12,2)],
           ASM = ASM[,seq(2,12,2)],
           colvec = myColor,
           plotSNP = FALSE) +
  labs(color = "Tissue")

ggdraw() +
  draw_plot(m2, x = 0, y = .6, width = .42, height = .4) +
  draw_plot(m1, x = .42, y = .6, width = .58, height = .4) +
  draw_plot(m3, x = 0, y = 0,width = 1, height = 0.6) +
  draw_plot_label(label = c("A", "B", "C"), size = 13,
                  x = c(0, 0.42, 0), y = c(1, 1, 0.6))

ggplot2::ggsave("curvesNscatters/MDSboth_newfunc.png", width = 10, height = 9)

#figure 6
source("custom_scoretracks.R")

#MEG3
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5)]
DAME <- GRanges(14, IRanges(101291540,101293480))
megcimp <- dame_track_forpap(DAME, 
                 window = 10, 
                 #positions = 400,
                 #derASM = derASM[,seq(2,12,2)],
                 ASM = ASM[,seq(2,12,2)],
                 colvec = myColor) +
  labs(color = "Tissue")
  

myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(3,8)]
megnon <- dame_track_forpap(DAME, 
                  window = 10, 
                  #positions = 400,
                  #derASM = derASM[,seq(2,12,2)],
                  ASM = ASM[,seq(1,11,2)],
                  colvec = myColor) +
  labs(color = "Tissue")

#h19
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5)]
DAME <- GRanges(11, IRanges(2021017,2021260))
h19cimp <- dame_track_forpap(DAME, 
                      window = 10, 
                      #positions = 400,
                      #derASM = derASM[,seq(2,12,2)],
                      ASM = ASM[,seq(2,12,2)],
                      colvec = myColor) +
  labs(color = "Tissue")

#gnas
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,5)]
DAME <- GRanges(20, IRanges(57425758,57428036))
gnascimp <- dame_track_forpap(DAME, 
                      #window = 10, 
                      positions = 700,
                      #derASM = derASM[,seq(2,12,2)],
                      ASM = ASM[,seq(2,12,2)],
                      colvec = myColor) +
  labs(color = "Tissue")

myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(3,8)]
gnasnon <- dame_track_forpap(DAME, 
                       #window = 10, 
                       positions = 700,
                       #derASM = derASM[,seq(2,12,2)],
                       ASM = ASM[,seq(1,11,2)],
                       colvec = myColor) +
  labs(color = "Tissue")

cowplot::plot_grid(megcimp, megnon, h19cimp, gnascimp, gnasnon, nrow = 5, 
                   ncol = 1, labels = c("A","","B", "C",""), align = "v")
ggplot2::ggsave(filename = "curvesNscatters/LOI_func.png",
                width = 10, height = 10)

#https://www.cell.com/action/showPdf?pii=S0167-7799%2818%2930115-X
#https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz125/5341422
