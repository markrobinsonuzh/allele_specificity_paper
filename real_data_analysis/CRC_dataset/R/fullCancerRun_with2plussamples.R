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
load("data/derASM_fullcancer2.RData")

#derASM <- GenomeInfoDb::sortSeqlevels(derASM) #only necessary for old calc_derivedasm()
#derASM <- sort(derASM)

#Filter
filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) >= 10
derASM <- derASM[filt,]

#Plot checks
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(2,3,5,8)]

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
mod <- model.matrix(~0+grp+samp)
mod <- mod[,-9] #because Coefficients not estimable: samp6_cimp
cont <- limma::makeContrasts(grpCRC_cimp-grpNORM_cimp, grpCRC_non-grpNORM_non, 
                             levels = mod)

#Get DAMEs
dames_cimp <- find_dames(derASM, mod, coef = 1, contrast = cont, maxGap = 100)
dames_noncimp <- find_dames(derASM, mod, coef = 2, contrast = cont, maxGap = 100)

#perms
dames_cimp <- find_dames(derASM, mod, coef = 1, contrast = cont, 
                            pvalAssign = "empirical", maxGap = 100) 

#use ASMstat
# refmeth <- assay(derASM, "ref.meth")
# altmeth <- assay(derASM, "alt.meth")
# 
# refcov <- assay(derASM, "ref.cov")
# altcov <- assay(derASM, "alt.cov")
# 
# prop <- (refmeth+altmeth)/(refcov+altcov)
# 
# ASMstat <- ((refmeth/refcov) - (altmeth/altcov)) /
#   sqrt(prop * (1 - prop) * ((1/refcov) + (1/altcov)))
# assay(derASM, "der.ASM") <- abs(ASMstat) #and run the above


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
colnames(derASM) <- c(paste0("C",1:6),paste0("N",1:6))
colnames(ASM) <- c(paste0("C",1:6),paste0("N",1:6))
m1 <- methyl_MDS_plot(derASM, group = metadata$V2, color = myColor, pointSize = 6, adj = 0.04) +
  theme(legend.position = "none")
m2 <- methyl_MDS_plot(ASM, group = metadata$V2, color = myColor, pointSize = 6, adj = 0.05) +
  theme(legend.position = "none")

#Paper figure
m4 <- cowplot::plot_grid(m2,m1, ncol=2, nrow = 2, labels = c("A","B","C"),
                         rel_widths = c(1,1), rel_heights = 1)

ggplot2::ggsave("curvesNscatters/MDSboth.png", m4, width = 8, height = 6)

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

#### build bigwigs ####

sapply(colnames(derASM), make_bigwig, scoreObj = derASM, folder =".." , 
       chromsizesFile = "/home/Shared_taupo/steph/reference/hg19.chrom.sizes.mod")

sapply(colnames(ASM), make_bigwig, scoreObj = ASM, folder =".." , 
       chromsizesFile = "/home/Shared_taupo/steph/reference/hg19.chrom.sizes.mod")
