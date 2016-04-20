# We use the sm_t matrix which contains the transformed ASM scores and the bump hunting methods to detect DAMEs.

# setwd("/home/Shared_taupo/data/seq/bisulphite_mirco/FASTQ/trimmed_t20l20/bump_hunting/full_bumphunting_real_data/with_bumphunter_package/")

# load sm_t
# load("real_sm_transformed.rda")

library(limma)
library(bumphunter)

############################################################

# The tuples in sm_t are already sorted by the median of the two CpG positions
# The rownames contain the chr, pos1, and pos2 per ASM score in the matrix
rows <- rownames(sm_t)
chr <- limma::strsplit2(rows, ".", fixed=T)[,1]
pos1 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,2])
pos2 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,3])
midpt <- floor((pos2 - pos1)/2)
pos <- pos1 + midpt

# moderated t-statistic
# mod is the design matrix. 0 is for normal, 1 is for adenoma
mod <- matrix(data=c(rep(1,13), c(0,0,1,1,1,1,1,1,1,1,1,0,1)), ncol = 2)
coeff <- 2 # the column in the design matrix to consider

fit <- limma::lmFit(sm_t, mod)

fit2 <- limma::eBayes(fit, proportion = 0.01)
beta <- fit2$t[, coeff]

# smoothing: by clusters on the genome called pns
# smoothing functions from bumphunter: loessByCluster and runmedByCluster

pns <- clusterMaker(chr, pos, maxGap = 300)

verbose <- TRUE
Q <- 0.9

smooth <- smoother(y = beta, x = pos, cluster = pns, smoothFunction = loessByCluster, 
                  verbose = verbose)
smoothed_beta <- smooth$fitted

# Detect DAMEs
dames <- regionFinder(x = smoothed_beta, chr = chr, pos = pos, cluster = pns, 
                      cutoff = quantile(abs(smoothed_beta),Q, na.rm=T), verbose = verbose)

save(dames, file="real_DAMEs.rda")

# number of regions DAMEfinder detects
# nrow(dames)
# [1] 22609

# Calculate a FDR using the null distribution of the premuted results (estimated FDR from permutaiton tests)
# We permute the labels of the samples on the smt_t matrix

# Get list of all possible permutations and remove the true one (We have 13 samples, 3 of which are normal_crypts)
combs <- combn(13, 3)
perm_list <- list()
for (i in 1:ncol(combs)) {
  v <- rep(1,13)
  v[combs[,i]] <- 0
  perm_list[[i]] <- v
}
w <- sapply(perm_list, function(x){all(x==c(0,0,1,1,1,1,1,1,1,1,1,0,1))}) 
w <- which(w==TRUE)
perm_list <- perm_list[-w]

# create a function that returns DAMEs given a permutation
pred_dames <- function(p) {
  # set design matrix
  mod <- matrix(data=c(rep(1,13),p), ncol = 2)
  # get t-stats
  fit <- limma::lmFit(sm_t, mod)
  fit2 <- limma::eBayes(fit, proportion = 0.01)
  beta <- fit2$t[, coeff]
  # do smoothing
  smooth <- smoother(y = beta, x = pos, cluster = pns, smoothFunction = loessByCluster, 
                    verbose = verbose)
  smoothed_beta <- smooth$fitted
  # predict regions
  dames <- regionFinder(x = smoothed_beta, chr = chr, pos = pos, cluster = pns, 
                        cutoff = quantile(abs(smoothed_beta),Q, na.rm=T), verbose = verbose)
  return(dames)
}

# Get observed DAMRs (true label)
obs_dames <- pred_dames(c(0,0,1,1,1,1,1,1,1,1,1,0,1))

# Get DAMEs with rest of permutations
perm_dames <- mclapply(perm_list, pred_dames, mc.cores=10)

# get values of the areas (measure of significance) of each region
perm_areas <- list()
for (i in 1:length(perm_dames)){
  perm_areas[[i]] <- perm_dames[[i]]$area
}

# calculate estimated p_value
# they are in order (lowest to highest) since the regions are soretd by area already
all_perm_areas <- unlist(perm_areas)
total_perm_areas <- length(all_perm_areas)
obs_dames$empirical_p_value <- rep(NA, nrow(obs_dames))
for (i in 1:nrow(obs_dames)) {
  w <- length(!is.na(which(obs_dames$area[i]>=all_perm_areas)))
  p <- (total_perm_areas-w)/total_perm_areas
  obs_dames[i,"empirical_p_value"] <- p
} 

save(obs_dames, file="real_DAMEs_with_fdr.rda")

# get list of DAMEs that have p-value <= 0.05
w <- which(obs_dames$empirical_p_value<=0.05)
obs_dames_fdr_0_0_5 <- obs_dames[w,]
save(obs_dames_fdr_0_0_5, file = "real_DAMEs_with_fdr_0_pt_0_5.rda")







