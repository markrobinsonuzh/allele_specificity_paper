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
library(ggplot2)

metadata <- read.table("data/fullCancerSamples.txt", stringsAsFactors = FALSE)
sample_names <- metadata$V1
#blacklist <- rtracklayer::import.bed("hg19-blacklist.v2.bed.gz")
#seqlevels(blacklist) <- gsub("chr", "", seqlevels(blacklist))

#### SNP mode ####

bam_files <- metadata$V3
vcf_files <- metadata$V4 
reference_file <- "/home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91/GRCh37.91.fa"

#Split reads and extract methylation according to allele
rds <- extract_bams(bam_files, vcf_files, sample_names, reference_file, cores=6, coverage = 2)
save(rds,file="data/extracted_bams.RData")
derASM <- calc_derivedasm(rds)

#save(derASM, file = "../data/derASM_fullCancer.RData")
load("data/derASM_fullcancer2.RData") #1,453,873

#derASM <- GenomeInfoDb::sortSeqlevels(derASM) #only necessary for old calc_derivedasm()
#derASM <- sort(derASM)

#Filter
filt <- c(rowSums(!is.na(assay(derASM, "der.ASM"))) >= 10
          & rowMeans(assay(derASM,"ref.cov") + assay(derASM,"alt.cov"), 
                     na.rm = TRUE) < 200)
derASM <- derASM[filt,] #55,717

#Plot checks
x <- assay(derASM,"der.ASM")
#x <- abs(ASMstat)

means <- rowMeans(x)
diffs <- apply(x, 1, function(w){mean(w[1:6]) - mean(w[7:12])})
var <- rowVars(x)
meancov <- rowMeans(assay(derASM,"ref.cov") + assay(derASM,"alt.cov"))
dd <- as.data.frame(cbind(var, means, diffs, meancov))

#MD plot
md <- ggplot(dd, aes(means, diffs)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  labs(x = "Average ASMsnp", y = "ASMsnp changes") +
  theme_bw()

#MV plot
MV <- ggplot(dd, aes(means, var)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  labs(x = "Average ASMsnp", y = "Variance ASMsnp") +
  theme_bw()

#var Vs cov
vacov <- ggplot(dd, aes(var, meancov)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  labs(x = "Variance ASMsnp", y = "Average coverage") +
  theme_bw()

#mean Vs cov
mecov <- ggplot(dd, aes(means, meancov)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  labs(x = "Average ASMsnp", y = "Average coverage") +
  theme_bw() 
#change axis names!
a <- cowplot::plot_grid(md, MV, vacov, mecov, nrow = 2, ncol = 2,
                        labels = c("A","B","C","D"))
ggplot2::ggsave(filename = "curvesNscatters/CRCdata_ASMsnp_diagnostics.png")
              

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
                            pvalAssign = "empirical", maxGap = 100, Q = 0.5)
dames_noncimp <- find_dames(derASM, mod, coef = 2, contrast = cont, 
                         pvalAssign = "empirical", maxGap = 100, Q = 0.5)

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
#    sqrt(prop * (1 - prop) * ((1/refcov) + (1/altcov)))



#### tuple mode ####

tuple_files <- metadata$V5
tuple_list <- read_tuples(files = tuple_files, sample_names, minCoverage = 5) 
ASM <- calc_asm(sampleList = tuple_list, transform = keepval) 
#keepval <- function(values){values * 1}
#use transform  = keepval for plotting
#save(ASM, file = "tupleASM_fullCancer_notrans.RData")
#save(ASM, file = "tupleASM_fullCancer.RData")
load("data/tupleASM_fullCancer.RData") #3,589,472

#Filter
filt <- c(rowSums(assay(ASM, "cov") >= 10 & !is.na(assay(ASM, "cov"))) >= 10
          & rowMeans(assay(ASM, "cov"), na.rm = TRUE) < 200)
ASM <- ASM[filt,] # 1,846,858

#MD and MV plots
x <- assay(ASM,"asm")
means <- rowMeans(x)
diffs <- apply(x, 1, function(w){mean(w[1:6]) - mean(w[7:12])})
var <- rowVars(x)
meancov <- rowMeans(assay(ASM,"cov"))
dd <- as.data.frame(cbind(var, means, diffs, meancov))

#MD plot
md <- ggplot(dd, aes(means, diffs)) + geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  labs(x = "Average ASMtuple", y = "ASMtuple changes") +
  theme_bw()

#MV plot
MV <- ggplot(dd, aes(means, var)) + geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  labs(x = "Average ASMtuple", y = "Variance ASMtuple") +
  theme_bw()

#var Vs cov
vacov <- ggplot(dd, aes(var, meancov)) +
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  labs(x = "Variance ASMtuple", y = "Average coverage") +
  theme_bw()

#mean Vs cov
mecov <- ggplot(dd, aes(means, meancov)) +
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  labs(x = "Average ASMtuple", y = "Average coverage") +
  theme_bw()

b <- cowplot::plot_grid(md, MV, vacov, mecov, nrow = 2, ncol = 2, 
                        labels = c("A","B","C","D"))
ggplot2::ggsave(filename = "curvesNscatters/CRCdata_ASMtuple_diagnostics.png", b)

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
dames_cimp <- find_dames(ASM, mod, contrast = cont, coef = 1, maxGap = 200) 
dames_noncimp <- find_dames(ASM, mod, contrast = cont,coef = 2, maxGap = 200) 

dames_cimp <- find_dames(ASM, mod, contrast = cont,coef = 1, maxGap = 200, pvalAssign = "empirical",Q = 0.5)
dames_noncimp <- find_dames(ASM, mod, contrast = cont,coef = 2, maxGap = 200, pvalAssign = "empirical", Q = 0.5)

#### build bigwigs ####

sapply(colnames(derASM), make_bigwig, scoreObj = derASM, folder =".." , 
       chromsizesFile = "/home/Shared_taupo/steph/reference/hg19.chrom.sizes.mod")

sapply(colnames(ASM), make_bigwig, scoreObj = ASM, folder =".." , 
       chromsizesFile = "/home/Shared_taupo/steph/reference/hg19.chrom.sizes.mod")
