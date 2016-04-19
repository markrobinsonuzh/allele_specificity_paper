# We calculate the ASM score for each tuple across the samples and transform the ASM scores with a square root transformation.

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

save(real_data, file="real_data.rda")

# get key of unique tuples
all_pos <- do.call("rbind", real_data)
all_pos <- all_pos[,c(1,3,4)]
key <- paste0(all_pos$chr,'.',all_pos$pos1, '.', all_pos$pos2)
key <- unique(key)

# length(key)
# [1] 2468506






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


