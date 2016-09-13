## Figure 6 -- plot TPR vs FDR for bsse and DAMEfinder
## Data: Simulated noisy data set

## rds files needed: "", "sim_ASM_score_matrix_noisy.rds", 
## extra file: 
## What we do: 

## Libraries
library(iCOBRA)
library(bsseq)
library(GenomicRanges)
library(bumphunter)


## A # Positions of the truely differentiated AMRs in the simulated data sat
dame_start <- c(93581084, 95234613, 100704363, 101004068, 101290524, 102226752)
dame_end <- c(93582797, 95240341, 100707067, 101005681, 101294433, 102228854)
dame_gr <- GRanges(seqnames = Rle(c("chr14"), c(length(dame_start))), 
                   ranges = IRanges(dame_start, end = dame_end), 
                   strand = rep("*", length(dame_start)))


## B # bsseq: do permutation tests and calculate empirical p-value of the DMRs predicted by bsseq



# permute the data and predict regions (necessary step to calculate the p_values) #

## Load the cov files (result of bismark)
cov_files <- dir(path = "../data/figure_6/noisy_cov_files",
                 pattern = "bismark.cov", full.names = TRUE, recursive = TRUE)
sample_names <- limma::strsplit2(cov_files, "_CpG_")[,1]
sample_names <- limma::strsplit2(sample_names, "-")[,1]
names(cov_files) <- sample_names

## Get bsseq object
samples <- read.bismark(cov_files, sample_names)

samples.fit <- BSmooth(samples, mc.cores=6, verbose=T) 
samples.cov <- getCoverage(samples.fit)

sample_conditions <- limma::strsplit2(cov_files, "/", fixed=T)[,8]

## Permute the sample labels and predict DMRs with bsseq
perm_list_dmr <- list(c("adenoma", "adenoma", "normals", "adenoma", "normals", "normals"),
                      c("adenoma", "adenoma", "normals", "normals", "adenoma", "normals"),
                      c("adenoma", "adenoma", "normals", "normals", "normals", "adenoma"),
                      c("adenoma", "normals", "adenoma", "adenoma", "normals", "normals"),
                      c("adenoma", "normals", "adenoma", "normals", "adenoma", "normals"),
                      c("adenoma", "normals", "adenoma", "normals", "normals", "adenoma"),
                      c("normals", "adenoma", "adenoma", "adenoma", "normals", "normals"),
                      c("normals", "adenoma", "adenoma", "normals", "adenoma", "normals"),
                      c("normals", "adenoma", "adenoma", "normals", "normals", "adenoma"))


## Function that predictss DMRs with bsseq
pred_dmrs <- function(p) {
  meta <- data.frame(sample=sample_names, condition=p)
  rownames(meta) <- sample_names
  # note: set at least 2 samples per condition
  keepLoci <- which(rowSums(samples.cov[, meta$condition=="adenoma"] >=4) >=2 & 
                      rowSums(samples.cov[, meta$condition=="normals", drop=F] >=4) >=2)
  
  samples.loci <- samples.fit[keepLoci,]
  
  samples.tstat <- BSmooth.tstat(samples.loci, group1=which(meta$condition=="adenoma"), 
                                 group2=which(meta$condition=="normals"), estimate.var="group2", 
                                 local.correct=T, verbose=T)
  
  # DMRs
  dmrs_quantil<-dmrFinder(samples.tstat, qcutoff=c(0.005,0.995))
  dmrs_subset<-subset(dmrs_quantil, n>=3 &abs(meanDiff) >= 0.1)
  return(dmrs_subset)
}

## Observed regions
obs_dmrs <- pred_dmrs(sample_conditions)

## Permuted regions
perm_dmrs <- mclapply(perm_list_dmr, pred_dmrs, mc.cores=9)

area_list <- list()
for (i in 1:length(perm_dmrs)){
  area_list[[i]] <- abs(perm_dmrs[[i]]$areaStat)
}

area_list[[length(area_list)+1]] <- abs(obs_dmrs$areaStat)
labels <- c(rep("perm", length(area_list)-1), "real")
names(area_list) <- labels

## Calculate empirical p-value
perm_areas <- unlist(area_list[1:(length(area_list)-1)])
obs_areas <- unlist(area_list[length(area_list)])
total_perm_areas <- length(perm_areas)
# add estimated p_value for each region
obs_dmrs$empirical_p_value <- rep(NA, nrow(obs_dmrs))
for (i in 1:length(obs_areas)) {
  w <- length(!is.na(which(obs_areas[i]>=perm_areas)))
  p <- (total_perm_areas-w)/total_perm_areas
  obs_dmrs[i,"empirical_p_value"] <- p
} 

## Set truth label
dmrs_gr <- makeGRangesFromDataFrame(obs_dmrs, keep.extra.columns = TRUE)
o <- countOverlaps(dmrs_gr, dame_gr, type = "any")
obs_dmrs$label <- o



