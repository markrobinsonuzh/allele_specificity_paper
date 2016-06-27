## Figure: 3 - ASM distribution plot
## Data: real data
## rda files needed: autosomal_df_ASM.rda, chrY_df_ASM.rda, chrX_df_ASM.rda, chrMT_df_ASM.rda
## What we do: Here we plot the distribution of the ASM score. We show the distribution for all samples for autosomal chr, 
## chrY, chrX, and chrMT

# load("autosomal_df_ASM.rda")
# load("chrY_df_ASM.rda")
# load("chrX_df_ASM.rda")
# load("chrMT_df_ASM.rda")

library(ggplot2)
library(grid)
library(gridExtra)

autosomal_df$sample <- as.factor(autosomal_df$sample)
autosomal_df$gender <- as.factor(autosomal_df$gender)
p_autosomal <- ggplot(autosomal_df, aes(x=sample, y=ASM_score, colour = gender, fill=gender)) + geom_violin() 

chrY_df$sample <- as.factor(chrY_df$sample)
chrY_df$gender <- as.factor(chrY_df$gender)
p_Y <- ggplot(chrY_df, aes(x=sample, y=ASM_score, colour = gender, fill=gender)) + geom_violin() 

chrX_df$sample <- as.factor(chrX_df$sample)
chrX_df$gender <- as.factor(chrX_df$gender)
p_X <- ggplot(chrX_df, aes(x=sample, y=ASM_score, colour = gender, fill=gender)) + geom_violin() 

chrMT_df$sample <- as.factor(chrMT_df$sample)
chrMT_df$gender <- as.factor(chrMT_df$gender)
p_MT <- ggplot(chrMT_df, aes(x=sample, y=ASM_score, colour = gender, fill=gender)) + geom_violin() 

pdf(file = "./ASM_violin_plots.pdf", w=10,h=8)
# par(mfrow=c(2,2))
grid.arrange(p_autosomal + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title="autosomal chr",x="Sample", y = "ASM Score") + coord_cartesian(ylim = c(-2, 4)), 
             p_Y + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title="chrY",x="Sample", y = "ASM Score") + coord_cartesian(ylim = c(-2, 4)), 
             p_X + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title="chrX",x="Sample", y = "ASM Score") + coord_cartesian(ylim = c(-2, 4)), 
             p_MT + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title="chrMT",x="Sample", y = "ASM Score") + coord_cartesian(ylim = c(-2, 4)),
             ncol = 2, nrow = 2)
dev.off()
