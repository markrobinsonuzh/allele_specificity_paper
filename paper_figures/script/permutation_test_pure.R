## Figure 5 -- plot distribution of DAMEs after doing permutation test
## Data: Simulated pure data set

## rds files needed: "sim_ASM_score_matrix_pure.rds", 
## extra file: 
## What we do: Permute the sample labels and use DAMEfinder to predict DAMEs. We compare to the DAMEs predicted on the true sample labels 

## Libraries
library(bumphunter)
library(limma)
library(parallel)

## Load data
sim_ASM_score_martix_pure<- readRDS("../data/sim_ASM_score_matrix_pure.rds")

## Remove positions where all adenoma (i.e. our simulated treatment) samples or all normal samples have NA values
w1 <- which(is.na(sim_ASM_score_martix_pure[,"adenoma1"])&is.na(sim_ASM_score_martix_pure[,"adenoma2"])&is.na(sim_ASM_score_martix_pure[,"adenoma3"]))
w2 <- which(is.na(sim_ASM_score_martix_pure[,"normal1"])&is.na(sim_ASM_score_martix_pure[,"normal2"])&is.na(sim_ASM_score_martix_pure[,"normal3"]))
sm <- sim_ASM_score_martix_pure[c(-w1,-w2),] 

## Order positions by median of each tuple
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

## Sqrt transfrom the data
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

## Parameters for bumphunting
coeff <- 2
Q <- 0.9
verbose <- TRUE

## Permutation List
perm_list <- list(c(1,1,0,1,0,0), c(1,1,0,0,1,0), c(1,1,0,0,0,1), 
                  c(1,0,1,1,0,0), c(1,0,1,0,1,0), c(1,0,1,0,0,1), 
                  c(0,1,1,1,0,0), c(0,1,1,0,1,0), c(0,1,1,0,0,1))


## Function that predicts DAMEs given a permutation
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
  odmrs <- regionFinder(x = smoothed_beta, chr = chr, pos = pos, cluster = pns, 
                        cutoff = quantile(abs(smoothed_beta),Q, na.rm=T), verbose = verbose)
  return(odmrs)
}

## Detect DAMEs on true labels
true_regs <- pred_regions(c(1,1,1,0,0,0))

## Detect DAMEs on permuted labels
perm_regs <- mclapply(perm_list, pred_regions, mc.cores=9)

## Get DAME areas
perm_area_list <- list()
for (i in 1:length(perm_regs)){
  perm_area_list[[i]] <- perm_regs[[i]]$area
}
true_area <- true_regs$area

## Calculate empirical p-values
perm_areas <- unlist(perm_area_list[1:(length(perm_area_list))])
total_perm_areas <- length(perm_areas)

true_regs$empirical_p_value <- rep(NA, nrow(true_regs))
for (i in 1:length(true_area)) {
  w <- length(!is.na(which(true_area[i]>=perm_areas)))
  p <- (total_perm_areas-w)/total_perm_areas
  true_regs[i,"empirical_p_value"] <- p
} 

## Calculate empirical FDRs per t (area). for t from 1 to 300. We had 9 permutations
fdr_df <- data.frame(t=seq(1:300), fdr=rep(NA,300))
for (i in 1:300) {
  # get number of discoveries >= the given t
  null_count <- length(which(perm_areas>=i))
  obs_count <- length(which(true_area>=i))
  fdr_df[i,"fdr"] <- (null_count/9)/obs_count
}

## At t=192 fdr<= 0.05
w <- which(fdr_df[,2]<=0.05)
head(fdr_df[w,])

## Plot boxplot of log area values in observed and permuted cases
a <- list(log(perm_areas, base = 10), log(true_area, base = 10))
b <- boxplot(a)
b$names <- c("perm", "real")
pdf("../figures/permutation_test_pure_boxplot.pdf", w=10,h=8)
bxp(b)
abline(a=log(192, base=10),b=0,col='red')
dev.off()




