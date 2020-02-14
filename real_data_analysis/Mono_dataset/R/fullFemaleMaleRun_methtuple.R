#########################################################################################
# Run tuple-mode on female-male files
#
# WGBS-seq data females Vs Males, BluePrint project
# NOTE: Files included here too big to include in repo
#
# Stephany Orjuela, January 2019
#########################################################################################

library(DAMEfinder)
library(SummarizedExperiment)
library(ggplot2)
library(tidyr)


DATA_PATH_DIR <- "/home/Shared_s3it/sorjuela/EGAdata/EGAD00001002523"

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

tuple_files <- sapply(c("C000S5/C000S5_full_qs.CG.2.tsv.gz",
                        "C001UY/C001UY_full_qs.CG.2.tsv.gz",
                        "S000RD/S000RD_full_qs.CG.2.tsv.gz",
                        "C0010K/C0010K_full_qs.CG.2.tsv.gz",
                        "C004SQ/C004SQ_full_qs.CG.2.tsv.gz",
                        "C005PS/C005PS_full_qs.CG.2.tsv.gz"), get_data_path)

sample_names <- c("C000S5", "C001UY", "S000RD", "C0010K", "C004SQ", "C005PS")

#Run DAMEfinder
tuple_list <- read_tuples(files = tuple_files, sample_names, minCoverage = 5)
ASM_mat <- calc_asm(sampleList = tuple_list) #6,802,057
#save(ASM_mat, file="malevsfem_ASM_mat.RData")
load("malevsfem_ASM_mat.RData")

#Filter
ASM_mat <- ASM_mat[rowSums(
  !is.na(assays(ASM_mat)[["cov"]]) &
    assays(ASM_mat)[["cov"]] >= 10) == BiocGenerics::ncol(ASM_mat),] #3,198,678

#Prepare for plot
get_e <- function(mat){
  
  ASM_x <- mat[seqnames(mat) == "X",]
  ASM_3 <- mat[seqnames(mat) == "3",]
  
  asmX <- as.data.frame(assays(ASM_x)[["asm"]])
  asm3 <- as.data.frame(assays(ASM_3)[["asm"]])
  colnames(asmX) <- colnames(asm3) <- c(1:6)
  
  eX <- tidyr::gather(asmX, sample, valueasm)
  eX$chr <- "chrX"
  
  e3 <- tidyr::gather(asm3, sample, valueasm)
  e3$chr <- "chr3"
  
  e <- rbind(eX,e3)
  e$Gender <- 0
  e$Gender <- ifelse(e$sample %in% 1:3, "Male", "Female")
  return(e)
}

e <- get_e(ASM_mat)

#Boxplots
myColor <- RColorBrewer::brewer.pal(9, "Set1")[3:4]
pfull <- ggplot(e) +
  geom_violin(aes(x=sample, y=valueasm, fill=Gender, color = Gender),
              trim = FALSE, adjust = 1.5)+ #, scale = "width") +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = "grey",
               alpha = 0,outlier.shape = NA) +
  scale_y_continuous(trans='sqrt') +
  scale_fill_manual(values = myColor) +
  scale_color_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15),
        legend.position = "none") +
  labs(y="ASMtuple",x="") +
  facet_grid(~chr)

#ggsave("curvesNscatters/test1.png",pfull)

#MDS
methyl_MDS_plot(ASM_mat, color = c(rep("Male",3), rep("Female",3)))

#beta Vs asm
# mm+0.5*um+0.5mu/total
# beta <- as.data.frame((assays(ASM_mat)[["MM"]]+
#                         (0.5*assays(ASM_mat)[["UM"]])+
#                         (0.5*assays(ASM_mat)[["MU"]])) / assays(ASM_mat)[["cov"]])
# ASM_meth <-ASM_mat
# assays(ASM_meth)[["asm"]] <- beta
# asm <- as.data.frame(assays(ASM_score_matrix)[["asm"]])
# asm <- as.data.frame(assays(ASM_score_matrix)[["cov"]])

#### Test for differences ####
males <- rowMeans(assay(ASM_x,"asm")[,1:3])
females <- rowMeans(assay(ASM_x,"asm")[,4:6])
var.test(females,males)

proms <- build_annotations(genome = 'hg38', annotations = "hg38_genes_promoters")

calc_diff <- function(ASM){
  cpgsites <- rowRanges(ASM)
  seqlevels(cpgsites) <- "chrX"
  over <- findOverlaps(cpgsites, proms)
  uniQ <- unique(queryHits(over))
  
  ASM.inproms <- ASM[uniQ,]
  males <- rowMeans(assay(ASM.inproms,"asm")[,1:3])
  females <- rowMeans(assay(ASM.inproms,"asm")[,4:6])
  #x <- t.test(females,males)
  x <- var.test(females,males)
  return(x)
}

calc_diff(ASM_x)

#### Boxplots only for promoters ####
library(annotatr)

proms <- build_annotations(genome = 'hg38', annotations = "hg38_genes_promoters")
cpgsites <- rowRanges(ASM_mat) # 3,198,678
seqlevels(cpgsites) <- paste0("chr", seqlevels(cpgsites))

over <- findOverlaps(cpgsites, proms) #1,666,247
uniQ <- unique(queryHits(over))

ASM_mat.inproms <- ASM_mat[uniQ,]

e <- get_e(ASM_mat.inproms)
e$chr <- ifelse(e$chr == "chrX", "chrX (promoters)", e$chr)
e$chr <- ifelse(e$chr == "chr3", "chr3 (promoters)", e$chr)

pprom <- ggplot(e) +
  geom_violin(aes(x=sample, y=valueasm, fill=Gender, color = Gender),
              trim = FALSE, adjust = 1.5) +# scale = "width") +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = "grey",
                alpha = 0,outlier.shape = NA) +
  scale_y_continuous(trans='sqrt') +
  scale_color_manual(values = myColor) +
  scale_fill_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15)) +
  labs(y="ASMtuple",x = "") +
  facet_grid(~chr)

#ggsave("curvesNscatters/test1.png",pprom)

#generate impr from other script, and use it here
p4 <- cowplot::plot_grid(pfull, pprom, impr, ncol=1, nrow = 3, 
                         labels = c("A","B","C"),
                         align = "v",
                         axis = "lr")
ggplot2::ggsave("curvesNscatters/chromboxplots_cow_noscalewidth.png", p4, 
                width = 10, height = 13)

#### make bigwigs ####

ASM_mat <- ASM_mat[seqnames(ASM_mat) %in% c(1:22, "X", "Y"),] #3,184,688

sapply(colnames(ASM_mat), make_bigwig, 
       scoreObj = ASM_mat, 
       folder = file.path(DATA_PATH_DIR, "bigwigs"), 
       chromsizesFile = "../../../Shared_taupo/steph/reference/hg38.chrom.sizes.mod")

