

#Get diffent scores with different weights

library(SummarizedExperiment)
library(DAMEfinder)
library(ggplot2)

metadata <- read.table("data/fullCancerSamples.txt", stringsAsFactors = FALSE)
sample_names <- metadata$V1
tuple_files <- metadata$V5

tuple_list <- read_tuples(files = tuple_files, sample_names, minCoverage = 5) 

#Filter
tuples_list_N1 <- tuple_list[7]
keep <- tuples_list_N1$NORM1_non$cov >= 10 & !is.na(tuples_list_N1$NORM1_non$cov)
tuples_list_N1$NORM1_non <- tuples_list_N1$NORM1_non[keep,]

#Set function to abs ASM_tuple
keepval <- function(values){abs(values * 1)}
myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(1:5,8:9)]

#Set gamma
betas <- c(0,0.5,0.8,1,5,10)

#Loop through all gammas
asms_b <- sapply(betas, function(i){
  ASM <- calc_asm(sampleList = tuples_list_N1, transform = keepval,
                beta = i, #penalize
                a = 0.2) #degree of departure
  asm_def <- as.vector(assay(ASM, "asm"))
  asm_def
})

colnames(asms_b) <- as.character(betas) 
asms_df <- as.data.frame(asms_b)
asms_df_long <- reshape2::melt(asms_df, value.name = "ASMtuple")
head(asms_df_long)

p1 <- ggplot(asms_df_long) + #test facets?
  geom_histogram(aes(ASMtuple, color = variable, fill = variable), alpha = 0.2) +
  scale_x_continuous(trans='sqrt') +
  scale_y_continuous(trans='sqrt') +
  facet_wrap(~variable) +
  labs(color = "gamma", fill = "gamma") +
  scale_color_manual(values = myColor) +
  scale_fill_manual(values = myColor) +
  theme_bw()
#ggsave("curvesNscatters/tupleASM_betas_hist.png", p1, width = 8, height = 6)

#set epsilon, cannot go above 0.5 since it's a value added or deleted from 0.5, 
#to result in a quantile 
epsil <- seq(0,0.5,0.1)

#Loop through all gammas
asms <- sapply(epsil, function(i){
  ASM <- calc_asm(sampleList = tuples_list_N1, transform = keepval,
                  beta = 0.5, #penalize
                  a = i) #degree of departure
  asm_def <- as.vector(assay(ASM, "asm"))
  asm_def
})

colnames(asms) <- as.character(epsil) 
asms_df <- as.data.frame(asms)
asms_df_long <- reshape2::melt(asms_df, value.name = "ASMtuple")
head(asms_df_long)

p2 <- ggplot(asms_df_long) + 
  geom_histogram(aes(ASMtuple, color = variable, fill = variable), alpha = 0.2) +
  scale_x_continuous(trans='sqrt') +
  scale_y_continuous(trans='sqrt') +
  facet_wrap(~variable) +
  labs(color = "epsilon", fill = "epsilon") +
  scale_color_manual(values = myColor) +
  scale_fill_manual(values = myColor) +
  theme_bw()

a <- cowplot::plot_grid(p1,p2, ncol=1, nrow = 2, labels = c("A","B"))
ggsave("curvesNscatters/tupleASM_weight_gammaANDepsilon.png", a, width = 6, height = 7.5)
