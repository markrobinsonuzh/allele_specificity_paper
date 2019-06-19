#########################################################################################
# Re-run full pipeline with non-cimp and cimp samples
#
# TBS-seq data CRCs Vs Norm
# NOTE: Files included here too big to include in repo
#
# Stephany Orjuela, January 2019
#########################################################################################

library(SummarizedExperiment)
library(DAMEfinder)

metadata <- read.table("../data/fullCancerSamples.txt", stringsAsFactors = FALSE)

#### SNP mode ####

bam_files <- metadata$V3
vcf_files <- metadata$V4 
sample_names <- metadata$V1
reference_file <- "/home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91/GRCh37.91.fa"

#Split reads and extract methylation according to allele
rds <- extract_bams(bam_files, vcf_files, sample_names, reference_file, cores=6)
derASM <- calc_derivedasm(rds)

#save(derASM, file = "../data/derASM_fullCancer.RData")
#load("../data/derASM_fullCancer.RData")

derASM <- GenomeInfoDb::sortSeqlevels(derASM) #only necessary for old calc_derivedasm()
derASM <- sort(derASM)

#Filter
filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) >= 10
derASM <- derASM[filt,]

#Plot checks
methyl_MDS_plot(derASM, color = metadata$V2)
methyl_MDS_plot(derASM, color = metadata$V2, top  = 500)
methyl_MDS_plot(derASM[,grp %in% c("NORM_cimp","NORM_non")], 
                color = grp[grp %in% c("NORM_cimp","NORM_non")], top = 10000)


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
#nice region
#9  99984229  99984408
dames_noncimp <- find_dames(derASM, mod, coef = 2, pvalAssign = "empirical", maxGap = 100) 

#Play with some regiones
which(start(derASM) == 99983700)
assay(derASM, "snp.table")[24208:24240,c(2,4,6,8,10,12)]
snp2 <- GRanges(9, IRanges(99983847, width = 1)) #este 
dame <- GRanges(9,IRanges(99983700,99983917)) #con este


#Plot methyl_circles 
#TODO: Change paths
methyl_circle_plot(snp2, vcfFile = gsub("/home/","/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/",metadata$V4[6]),
                   bamFile = gsub("/home/","/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/",metadata$V3[6]), 
                   refFile = gsub("/home/","/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/",reference_file),
                   sampleName = "CRC6",
                   dame = dame,
                   pointSize = 2)
                   
methyl_circle_plot(snp2, vcfFile = gsub("/home/","/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/",metadata$V4[12]),
                   bamFile = gsub("/home/","/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/",metadata$V3[12]), 
                   refFile = gsub("/home/","/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/",reference_file),
                   sampleName = "NORM6",
                   pointSize = 2)



#### tuple mode ####

tuple_files <- metadata$V5
tuple_list <- read_tuples(files = tuple_files, sample_names, min_coverage = 5)
ASM <- calc_asm(sample_list = tuple_list)
#save(ASM, file = "tupleASM_fullCancer.RData")
#load("../tupleASM_fullCancer.RData")

#Filter
filt <- c(rowSums(assay(ASM, "cov") >= 10 & !is.na(assay(ASM, "cov"))) >= 10)
ASM <- ASM[filt,] #2,015,001, 1,576,545

#MDS
methyl_MDS_plot(ASM, color = metadata$V2)

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
mod <- model.matrix(~grp+samp)
mod <- mod[,-9]

dames_noncimp <- find_dames(ASM, mod, coef = 2, maxGap = 200)

#### build bigwigs ####

sapply(colnames(derASM), make_bigwig, score.obj = derASM, folder =".." , 
       chromsizes.file = "/home/Shared_taupo/steph/reference/hg19.chrom.sizes.mod")

sapply(colnames(ASM), make_bigwig, score.obj = ASM, folder =".." , 
       chromsizes.file = "/home/Shared_taupo/steph/reference/hg19.chrom.sizes.mod")


