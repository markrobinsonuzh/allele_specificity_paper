#########################################################################################
# Re-run full pipeline with non-cimp and cimp samples, including ASMsnp and ASMtuple
#
# TBS-seq data CRCs Vs Norm
# NOTE: Files included here too big to include in repo
#
# Stephany Orjuela, January 2019
#########################################################################################

library(SummarizedExperiment)
library(DAMEfinder)

metadata <- read.table("data/fullCancerSamples.txt", stringsAsFactors = FALSE)

#### SNP mode ####

bam_files <- metadata$V3
vcf_files <- metadata$V4 
sample_names <- metadata$V1
reference_file <- "/home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91/GRCh37.91.fa"

#Split reads and extract methylation according to allele
rds <- extract_bams(bam_files, vcf_files, sample_names, reference_file, cores=6, coverage = 2)
save(rds,file="data/extracted_bams.RData")
derASM <- calc_derivedasm(rds)

#save(derASM, file = "../data/derASM_fullCancer.RData")
load("data/derASM_fullCancer.RData")

derASM <- GenomeInfoDb::sortSeqlevels(derASM) #only necessary for old calc_derivedasm()
derASM <- sort(derASM)

#Filter
filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) >= 10
derASM <- derASM[filt,]

#Plot checks
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,3,5,8)]
m1 <- methyl_MDS_plot(derASM, group = metadata$V2, color = myColor)

x <- assay(derASM,"der.ASM")
#xsub <- x[,grp %in% c("CRC_cimp","NORM_cimp")]
means <- rowMeans(x)
diffs <- apply(x, 1, function(w){mean(w[1:6]) - mean(w[7:12])})
var <- rowVars(x)
dd <- as.data.frame(cbind(var, means, diffs))

#MD plot
ggplot(dd, aes(means, diffs)) + geom_point(alpha = 0.2) +
  theme_bw()#+ ylim(c(-0.4,0.4))

#MV plot
ggplot(dd, aes(means, var)) + geom_point(alpha = 0.2) + theme_bw()

#Set design
grp <- factor(metadata$V2)
grp <- relevel(grp, "NORM_cimp")
samp <- gsub("CRC|NORM","", metadata$V1)
mod <- model.matrix(~grp+samp)
mod <- mod[,-9] #because Coefficients not estimable: samp6_cimp

#Get DAMEs
dames_noncimp <- find_dames(derASM, mod, coef = 2, maxGap = 100)
dames_noncimp <- find_dames(derASM, mod, coef = 2, pvalAssign = "empirical", maxGap = 100) 

#Plot methyl_circles 
snp2 <- GRanges(9, IRanges(99983847, width = 1)) #este 
dame <- GRanges(9,IRanges(99983700,99983917)) #con este

path <- "/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/"

m1 <- methyl_circle_plot(snp2, vcfFile = gsub("/home/",path, metadata$V4[2]),
                   bamFile = gsub("/home/",path,metadata$V3[2]), 
                   refFile = gsub("/home/",path,reference_file),
                   sampleName = "CRC2",
                   dame = dame,
                   pointSize = 2)
                   
m2 <- methyl_circle_plot(snp2, vcfFile = gsub("/home/",path,metadata$V4[8]),
                   bamFile = gsub("/home/",path,metadata$V3[8]), 
                   refFile = gsub("/home/",path,reference_file),
                   sampleName = "NORM2",
                   dame = dame,
                   pointSize = 2)

m4 <- cowplot::plot_grid(m2,m1, ncol=1, nrow = 2, labels = c("A","B"))
ggplot2::ggsave("MCircle_plots/MethylcirclesSNP.png", m4, width = 12, height = 10)



#### tuple mode ####

tuple_files <- metadata$V5
tuple_list <- read_tuples(files = tuple_files, sample_names, min_coverage = 5)
ASM <- calc_asm(sample_list = tuple_list)
#save(ASM, file = "tupleASM_fullCancer.RData")
load("data/tupleASM_fullCancer.RData")

#Filter
filt <- c(rowSums(assay(ASM, "cov") >= 10 & !is.na(assay(ASM, "cov"))) >= 10)
ASM <- ASM[filt,] #2,015,001, 1,849,831

#MDS
m2 <- methyl_MDS_plot(ASM, group = metadata$V2, color = myColor, adj = 0.03)

#full paper figure TODO: fix for new MDS function
m1 <- m1 + theme(text = element_text(size = 15))
m2 <- m2 +  theme(legend.position = "none", text = element_text(size = 15))
m4 <- cowplot::plot_grid(m2,m1, ncol=2, nrow = 2, labels = c("A","B","C"),
                         rel_widths = c(1,1.3), rel_heights = 1)
ggplot2::ggsave("curvesNscatters/MDSboth.png", m4, width = 12, height = 10)

#MD and MV plots
x <- assay(ASM,"asm")
means <- rowMeans(x)
var <- rowVars(x)
diffs <- apply(x, 1, function(w){mean(w[1:6]) - mean(w[7:12])})
dd <- as.data.frame(cbind(diffs, means))
ggplot(dd, aes(means, diffs)) + geom_point(alpha=0.2) + theme_bw()

#Set design
grp <- factor(metadata$V2)
grp <- relevel(grp, "NORM_cimp")
samp <- gsub("CRC|NORM","", metadata$V1)
mod <- model.matrix(~0+grp+samp)
mod <- mod[,-9]

#set contrast
cont <- limma::makeContrasts(grpCRC_cimp-grpNORM_cimp, grpCRC_non-grpNORM_non, 
                             levels = mod)

#get DAMEs
dames_cimp <- find_dames(ASM, mod, contrast = cont, coef = 1, maxGap = 200) #4037
dames_noncimp <- find_dames(ASM, mod, contrast = cont,coef = 2, maxGap = 200) #260
dames_noncimp <- dames_noncimp[dames_noncimp$FDR <= 0.05 & dames_noncimp$clusterL > 1,] #231

dames_cimp <- find_dames(ASM, mod, contrast = cont,coef = 1, maxGap = 200, pvalAssign = "empirical") #0
dames_noncimp <- find_dames(ASM, mod, contrast = cont,coef = 2, maxGap = 200, pvalAssign = "empirical")#0

#Genes with known LOI
MEG3 <- GRanges(14, IRanges(101245747,101327368)) 
GNAS <- GRanges(20, IRanges(57414773,57486247))
H19 <- GRanges(11,IRanges(2016406,2022700))
IGF2 <- GRanges(11,IRanges(2150342,2162468))


#Some methylcircles
trimmedMEG3 <- GRanges(12,IRanges(98850698,98851011))   
m1 <- methyl_circle_plotCpG(trimmedMEG3,
                      bamFile = gsub("/home/",path,metadata$V3[2]), 
                      refFile = gsub("/home/",path,reference_file),
                      dame = trimmedMEG3,
                      pointSize = 1,
                      order = TRUE)
m2 <- methyl_circle_plotCpG(trimmedMEG3,
                            bamFile = gsub("/home/",path,metadata$V3[8]), 
                            refFile = gsub("/home/",path,reference_file),
                            dame = trimmedMEG3,
                            pointSize = 1,
                            order = TRUE)

m4 <- cowplot::plot_grid(m2,m1, ncol=1, nrow = 2, labels = c("A","B"))
ggplot2::ggsave("MCircle_plots/methylCirclesTopTupleDAMECIMP.png", m4, width = 10, 
                height = 12)

#### build bigwigs ####

sapply(colnames(derASM), make_bigwig, score.obj = derASM, folder =".." , 
       chromsizes.file = "/home/Shared_taupo/steph/reference/hg19.chrom.sizes.mod")

sapply(colnames(ASM), make_bigwig, score.obj = ASM, folder =".." , 
       chromsizes.file = "/home/Shared_taupo/steph/reference/hg19.chrom.sizes.mod")


