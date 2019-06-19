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

#Since the files are too big, I only run the chromosomes I want to plot (see BASH folder)

#chr3
tuple_files <- sapply(c("C000S5/C000S5_single_3subset_qs.CG.2.tsv.gz",
                      "C001UY/C001UY_single_3subset_qs.CG.2.tsv.gz",
                      "S000RD/S000RD_single_3subset_qs.CG.2.tsv.gz",
                      "C0010K/C0010K_single_3subset_qs.CG.2.tsv.gz",
                      "C004SQ/C004SQ_single_3subset_qs.CG.2.tsv.gz",
                      "C005PS/C005PS_single_3subset_qs.CG.2.tsv.gz"), get_data_path)

#chrX
tuple_files <- sapply(c("C000S5/C000S5_Xsubset_qs.CG.2.tsv.gz",
                        "C001UY/C001UY_Xsubset_qs.CG.2.tsv.gz",
                        #"S000RD/S000RD_Xsubset_qs.CG.2.tsv.gz",
                        "S000RD/S000RD_single_Xsubset_qs.CG.2.tsv.gz",
                        "C0010K/C0010K_Xsubset_qs.CG.2.tsv.gz",
                        "C004SQ/C004SQ_Xsubset_qs.CG.2.tsv.gz",
                        "C005PS/C005PS_Xsubset_qs.CG.2.tsv.gz"), get_data_path)

sample_names <- c("C000S5", "C001UY", "S000RD", "C0010K", "C004SQ", "C005PS")

#Run DAMEfinder for each chromosome
tuple_list <- read_tuples(files = tuple_files, sample_names, min_coverage = 5)
ASM_mat <- calc_asm(sample_list = tuple_list)
#save(ASM_mat, "malevsfem_chrom3_ASM_mat.RData")


#Filter
#load("malevsfem_chromX_ASM_mat_singleTest.RData") #317,179
ASM_x <- ASM_mat[rowSums(
  !is.na(assays(ASM_mat)[["cov"]]) &
    assays(ASM_mat)[["cov"]] >= 10) == BiocGenerics::ncol(ASM_mat),] #X:18,365

#load("malevsfem_chrom3_ASM_mat.RData") #352,259
ASM_mat <- ASM_mat[rowSums(
  !is.na(assays(ASM_mat)[["cov"]]) &
    assays(ASM_mat)[["cov"]] >= 10) == BiocGenerics::ncol(ASM_mat),] #3:166,978

ASM_mat <- rbind(ASM_mat,ASM_x) #185,343


#Prepare for plot
asm <- as.data.frame(assays(ASM_mat)[["asm"]])
colnames(asm) <- c(1:6)

e <- tidyr::gather(asm, sample, valueasm)
e$Gender <- 0
e$Gender <- ifelse(e$sample %in% 1:3, "Male", "Female")
e$chr <- 0 

for(i in 1:6){
  e$chr[e$sample == i][1:166978] <- "chr3"
  e$chr[e$sample == i][166979:185343] <- "chrX"
}

#Boxplots
myColor <- RColorBrewer::brewer.pal(9, "Set1")[3:4]
greys <- RColorBrewer::brewer.pal(9, "Greys")[7]
ggplot(e) +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = greys) +
  scale_fill_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15)) +
  ylab("ASMtuple") +
  facet_grid(~chr)


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

#asm <- abs(asm)
#w <- !is.na(rowSums(beta))
#beta <- beta[w,]

# myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
# #olorRampPalette(bre[4:9])(15)
# 
# #myColor_scale_fill <- scale_fill_gradientn(colours = myColor)
# 
# ggplot(de, aes(x=value, y=valueasm)) +
#   #stat_binhex(bins = 50) +
#   #myColor_scale_fill +
#   geom_point(alpha = 1/10) +
#   theme_light() +
#   #geom_density2d(colour = "black") +
#   facet_wrap(~sample)


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
cpgsites <- rowRanges(ASM_mat)
seqlevels(cpgsites) <- c("chr3","chrX")

over <- findOverlaps(cpgsites, proms)
uniQ <- unique(queryHits(over))

ASM_mat.inproms <- ASM_mat[uniQ,]
asm <- as.data.frame(assays(ASM_mat.inproms)[["asm"]])
colnames(asm) <- c(1:6)

e <- tidyr::gather(asm, sample, valueasm)
e$Gender <- 0
e$Gender <- ifelse(e$sample %in% 1:3, "Male", "Female")
e$chr <- 0 

#To check:
seqnames(rowRanges(ASM_mat.inproms))

for(i in 1:6){
  e$chr[e$sample == i][1:27613] <- "chr3"
  e$chr[e$sample == i][27614:29349] <- "chrX"
}

ggplot(e) +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = greys) +
  scale_fill_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15)) +
  ylab("ASMtuple") +
  facet_grid(~chr)

