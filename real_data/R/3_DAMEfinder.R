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

# Calculate a FDR using the null distribution of the premuted results (estimated FDR from permutaiton tests)
# We permute the labels of the samples on the smt_t matrix

# Get list of all possible permutations and remove the true one





