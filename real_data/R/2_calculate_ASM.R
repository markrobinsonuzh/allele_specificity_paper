# We calculate the ASM score for each tuple across the samples and transform the ASM scores with a square root transformation. The final matrix is called the sm_t matrix. It has the tuples sorted by position of the median per tuple and contains the transformed ASM scores.

library(parallel)
library(data.table)
library(LearnBayes)

# setwd("/home/Shared/data/seq/bisulphite_mirco/FASTQ/trimmed_t20l20/bump_hunting/full_bumphunting_real_data/with_bumphunter_package")

########################################################################

# load("real_methtuple_list_o_final.rda")

# add 1 to every count and calculate the log odds ratio per tuple
real_data <- mclapply(real_methtuple_list_o_final, process, mc.cores = length(real_methtuple_list_o_final))

# set parameters for ASM score. 'beta' values for beta binomial model and 'a' the distance from 0.5 we consider the probability of
beta <- 0.5
a <- 0.2

# calculate ASM score per tuple per sample
for (i in 1:length(real_data)) {
  df <- real_data[[i]]
  df$ASM_score <- calcScore(beta,a,df)
  real_data[[i]] <- df
}

# save(real_data, file="real_data.rda")

# get key of unique tuples
all_pos <- do.call("rbind", real_data)
all_pos <- all_pos[,c(1,3,4)]
key <- paste0(all_pos$chr,'.',all_pos$pos1, '.', all_pos$pos2)
key <- unique(key)

# length(key)
# [1] 2468506

# get matrix of ASM scores across all samples
real_score_matrix <- matrix(data=NA, nrow=length(key), ncol=length(real_data))
colnames(real_score_matrix) <- paste0(names(real_data), c("-normal", "-normal", "-adenoma", 
                                                          "-adenoma", "-adenoma", "-adenoma", 
                                                          "-adenoma", "-adenoma", "-adenoma", 
                                                          "-adenoma", "-adenoma", "-normal",
                                                          "-adenoma"))
rownames(real_score_matrix) <- key
for (i in 1:length(real_data)) {
  df <- real_data[[i]]
  key_s <- paste0(df$chr, '.', df$pos1, '.', df$pos2)
  m <- match(key, key_s)
  ind <- m[which(!is.na(m))]
  ind_m <- which(!is.na(m))
  real_score_matrix[ind_m,i] <- df$ASM_score[ind]
}

# save(real_score_matrix, file="real_score_matrix.rda")

# nrow(real_score_matrix)
# [1] 2468506

# remove the positions where all the normal samples have NA values or all the adenoma samples have NA values and store as sm matrix

n <- grep("normal", colnames(real_score_matrix))
a <- grep("adenoma", colnames(real_score_matrix))

norm_scores <- real_score_matrix[,n]
w_n <- which(rowSums(is.na(norm_scores))==3)

adenoma_scores <- real_score_matrix[,a]
w_a <- which(rowSums(is.na(adenoma_scores))==10)

sm <- real_score_matrix[c(-w_n,-w_a),] 

# nrow(sm)
# [1] 1638587

# set median of each tuple as the genomic position
pos <- rownames(sm)
chr <- limma::strsplit2(pos, ".", fixed=T)[,1]
pos1 <- as.numeric(limma::strsplit2(pos, ".", fixed=T)[,2])
pos2 <- as.numeric(limma::strsplit2(pos, ".", fixed=T)[,3])
midpt <- floor((pos2 - pos1)/2)
median_pos <- pos1 + midpt

# sort the score matrix by chr and median position (important for regionFinder and bumphunting)
pos_df <- data.frame(chr=chr, pos=median_pos)
o <- order(pos_df[,"chr"], pos_df[,"pos"])
sm <- sm[o,]
median_pos <- median_pos[o]
chr <- chr[o]

# save the sm matrix
save(sm, file="real_sm.rda")

# transform the ASM scores. We do this to have a more stable mean-variance relationship when we use limma in a later step.
sm_t <- matrix(data=NA, nrow=nrow(sm), ncol=ncol(sm))
colnames(sm_t) <- colnames(sm)
rownames(sm_t) <- rownames(sm)
for (i in 1:ncol(sm_t)) {
  sm_t[,i] <- modulus_sqrt(sm[,i])
}

save(sm_t, file="real_sm_transformed.rda")






## list of functions

# adds one to every cell and calculate log odds ratio per tuple in a given data frame
process <- function(s) {
  s$strand <- "*"
  s$MM <- s$MM + 1
  s$MU <- s$MU + 1
  s$UM <- s$UM + 1
  s$UU <- s$UU + 1
  s$cov <- s$cov + 4
  # calc log odds ratio
  ratio <- (s$MM*s$UU)/(s$MU*s$UM)
  s$log_odds_ratio <-  log(ratio, base=10)
  return(s)
}


# calculate the weight per site given beta and a
calcWeight <- function(beta, MM, UU, a) {
  
  # set quantile range
  theta = seq(0.005, 0.995, length = 1000)
  
  # calculate pbeta (note that this is the cumulative prob distr). To get the densities, use dbeta.
  posterior <- pbeta(theta, beta+MM, beta+UU, lower.tail = T)
  
  # get prob(0.3<b<0.7, ie a=0.2): our new weight
  w_1 <- which(theta<(0.5-a))
  index_1 <- w_1[length(w_1)] + 1
  
  w_2 <- which(theta<(0.5+a))
  index_2 <- w_2[length(w_2)]
  
  weight <- posterior[index_2] - posterior[index_1]
  
  return(weight)
  
}


# get a vector of the posterior scores for all positions given beta and a
calcScore <- function(beta, a, df) {
  
  weights <- mcmapply(calcWeight, 1, df$MM, df$UU, a, mc.cores = 7)
  score_vector <- df$log_odds_ratio*weights
  return(score_vector)
  
}

# modulus sqrt function to  transform the data 
modulus_sqrt <- function(values) {
  t_values <- abs(values)
  ret <- sign(values)*sqrt(t_values)
  return(ret)
}


