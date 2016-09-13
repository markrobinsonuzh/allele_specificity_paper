## Figure 6 -- plot TPR vs FDR for bsse and DAMEfinder
## Data: Simulated noisy data set

## What we do: Look at performance of bsseq vs DAMEfinder on simulated data. We calculate the empirical p-value for each regions the two
## methods predict using permutation tests.

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


## C # DAMEfinder: do permutation tests and calculate empirical p-value of the DAMEs predicted by DAMEfinder

## Load ASM score matrix
sim_score_matrix <- readRDS("../data/sim_ASM_score_matrix_noisy.rds")

## Remove CpG positions where all normals or all adenoma are NA
w1 <- which(is.na(sim_score_matrix[,"adenoma1"])&is.na(sim_score_matrix[,"adenoma2"])&is.na(sim_score_matrix[,"adenoma3"]))
w2 <- which(is.na(sim_score_matrix[,"normal1"])&is.na(sim_score_matrix[,"normal2"])&is.na(sim_score_matrix[,"normal3"]))
sm <- sim_score_matrix[c(-w1,-w2),] 

## Order by position
rows <- rownames(sm)
chr <- limma::strsplit2(rows, ".", fixed=T)[,1]
pos1 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,2])
pos2 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,3])
midpt <- floor((pos2 - pos1)/2)
pos <- pos1 + midpt

pos_df <- data.frame(chr=chr, pos=pos)
o <- order(pos_df[,"chr"], pos_df[,"pos"])
sm <- sm[o,]
pos <- pos[o]
chr <- chr[o]

## Generate clusters
pns <- clusterMaker(chr, pos, maxGap = 300)

# Sqrt transfrom the data
modulus_sqrt <- function(values) {
  t_values <- abs(values)
  ret <- sign(values)*sqrt(t_values)
  return(ret)
}

sm_t <- matrix(data=NA, nrow=nrow(sm), ncol=ncol(sm))
colnames(sm_t) <- colnames(sm)
rownames(sm_t) <- rownames(sm)
for (i in 1:ncol(sm_t)) {
  sm_t[,i] <- modulus_sqrt(sm[,i])
}

coeff <- 2
Q <- 0.9
verbose <- TRUE

## Permutation test to calculate empirical p-value
perm_list <- list(c(1,1,0,1,0,0), c(1,1,0,0,1,0), c(1,1,0,0,0,1), 
                  c(1,0,1,1,0,0), c(1,0,1,0,1,0), c(1,0,1,0,0,1), 
                  c(0,1,1,1,0,0), c(0,1,1,0,1,0), c(0,1,1,0,0,1))

## Function that predicts DAMEs given a permutation vector
pred_regions <- function(p) {
  # set design matrix
  mod <- matrix(data=c(1,1,1,1,1,1,p), ncol = 2)
  # get t-stats
  fit <- limma::lmFit(sm_t, mod)
  fit2 <- limma::eBayes(fit, proportion = 0.01, robust = TRUE)
  wald <- fit2$t[, coeff]
  # do smoothing
  swald <- smoother(y = wald, x = pos, cluster = pns, smoothFunction = loessByCluster, 
                    verbose = verbose)
  smoothed_beta <- swald$fitted
  # predict regions
  dames <- regionFinder(x = smoothed_beta, chr = chr, pos = pos, cluster = pns, 
                        cutoff = quantile(abs(smoothed_beta),Q, na.rm=T), verbose = verbose)
  return(dames)
}

## Observed
obs_regs <- pred_regions(c(1,1,1,0,0,0))

## Permuted
perm_regs <- mclapply(perm_list, pred_regions, mc.cores=9)

## Areas and empirical p-values
area_list <- list()
for (i in 1:length(perm_regs)){
  area_list[[i]] <- perm_regs[[i]]$area
}
area_list[[length(area_list)+1]] <- obs_regs$area
labels <- c(rep("perm", length(area_list)-1), "real")
names(area_list) <- labels

perm_areas <- unlist(area_list[1:(length(area_list)-1)])
obs_areas <- unlist(area_list[length(area_list)])
total_perm_areas <- length(perm_areas)
obs_regs$empirical_p_value <- rep(NA, nrow(obs_regs))
for (i in 1:length(obs_areas)) {
  w <- length(!is.na(which(obs_areas[i]>=perm_areas)))
  p <- (total_perm_areas-w)/total_perm_areas
  obs_regs[i,"empirical_p_value"] <- p
} 

## Set truth label
obs_gr <- makeGRangesFromDataFrame(obs_regs, keep.extra.columns = TRUE)
o <- countOverlaps(obs_gr, dame_gr, type = "any")
obs_regs$label <- o


## D # iCOBRA

bsseq_dmrs <- paste0(obs_dmrs$chr,".",obs_dmrs$start,".",obs_dmrs$end)
DAMEfinder_dames <- paste0(obs_regs$chr,".",obs_regs$start,".",obs_regs$end)

p_values <- data.frame(p_value_bsseq = c(obs_dmrs$empirical_p_value, rep(NA,length(DAMEfinder_dames))), 
                       p_value_DAMEfinder = c(rep(NA, length(bsseq_dmrs)),obs_regs$empirical_p_value))
colnames(p_values) <- c("bsseq", "DAMEfinder")
rownames(p_values) <- c(bsseq_dmrs, DAMEfinder_dames)


truth <- data.frame(truth=c(obs_dmrs$label,obs_regs$label))
colnames(truth) <- "status"
rownames(truth) <- c(bsseq_dmrs, DAMEfinder_dames)

cobra_data <- COBRAData(pval = p_values, truth=truth)
cobra_data <- calculate_adjp(cobra_data)

cobraperf <- calculate_performance(cobra_data, binary_truth = "status", onlyshared=TRUE)

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", facetted = TRUE)

pdf("../figures/figure_6/tpr_fdr_bsseq_DAMEfinder_noisy.pdf", w=10, h=8)
plot_fdrtprcurve(cobraplot)
dev.off()

