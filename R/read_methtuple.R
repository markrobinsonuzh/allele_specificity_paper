# We read in the methtuple files of the real data set and do some filtering on the tuples we want to keep. We set a minimum coverage
# of 10 per tuple and a maxGap of 150 base pairs between two CpG sites in a tuple. We also remove positions that correspond to SNPs.
# The SNP database for hg19 was used to obtain a list of SNPs. "all_snps_hg19.rda" is an object containing a list of the SNPs (chr and position).

library("data.table")
library("parallel")

# setwd("/home/Shared_taupo/data/seq/bisulphite_mirco/FASTQ/trimmed_t20l20/bump_hunting/full_bumphunting_real_data/with_bumphunter_package")

#########################################################

min_cov <- 10
max_gap <- 150

# Read in the methtuple files as data frames into a list called real_methtuple_list
dir <- "/home/Shared_taupo/data/seq/bisulphite_mirco/FASTQ/trimmed_t20l20"
real_methtuple_files <- dir(path = dir, pattern = "sorted.CG.2.tsv.gz", full.names = TRUE, recursive = TRUE)
nm <- limma::strsplit2(real_methtuple_files, "QS_NGS")[,1]
nm <- limma::strsplit2(nm, "20/")[,2]
nm <- limma::strsplit2(nm, "/")[,1]
names(real_methtuple_files) <- nm

real_methtuple_list <- mclapply(real_methtuple_files, function(file){read.table(file, header=T)}, mc.cores=length(real_methtuple_files))

# order each data frame by chr, pos1 and pos2, add cov column, and add tuple distance column
real_methtuple_list_o <- list()
for (i in 1:length(real_methtuple_list)) {
  df <- real_methtuple_list[[i]]
  df$cov <- df$MM + df$UU + df$UM + df$MU
  df$inter_dist <- df$pos2 - df$pos1
  df <- df[order(df[,1], df[,3], df[,4]),]
  real_methtuple_list_o[[i]] <- df
  names(real_methtuple_list_o)[i] <- names(real_methtuple_list)[i]
}

# filter out by coverage and maxGap
real_methtuple_list_o_f <- list()
for (i in 1:length(real_methtuple_list_o)) {
  df <- real_methtuple_list_o[[i]]
  w <- which(df$cov >= min_cov & df$inter_dist <= max_gap)
  real_methtuple_list_o_f[[i]] <- df[w,]
  names(real_methtuple_list_o_f)[i] <- names(real_methtuple_list_o)[i]
}

# remove SNPs
load("all_snps_hg19.rda")
snp_key <- paste0("chr", all_snps$chr,".", all_snps$pos)

remove_snps <- function(df) {
    key_1 <- paste0(df$chr,".", df$pos1)
    key_2 <- paste0(df$chr,".", df$pos2)
    
    m1 <- match (snp_key, key_1)
    w1 <- m1[!is.na(m1)]
    
    m2 <- match (snp_key, key_2)
    w2 <- m2[!is.na(m2)]
    
    return(df[-c(w1,w2),])
}

real_methtuple_list_o_final <- mclapply(real_methtuple_list_o_f, remove_snps, mc.cores=length(real_methtuple_list_o_f))

# save(real_methtuple_list_o_final, file="real_methtuple_list_o_final.rda")




