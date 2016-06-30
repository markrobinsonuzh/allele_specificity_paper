## Figure: 3 - ASM distribution plot
## Data: real data
## rda files needed: autosomal_df_ASM.rda, chrY_df_ASM.rda, chrX_df_ASM.rda, chrMT_df_ASM.rda
## What we do: Here we plot the distribution of the ASM score. We show the distribution for all samples for autosomal chr, 
## chrY, chrX, and chrMT

load("../data/autosomal_df_ASM.rda")
load("../data/chrY_df_ASM.rda")
load("../data/chrX_df_ASM.rda")
#load("chrMT_df_ASM.rda")

library(ggplot2)

df <- rbind( cbind(autosomal_df, chrom="autosomal"),
             cbind(chrX_df, chrom="X"),
             cbind(chrY_df, chrom="Y"))

pdf("../figures/ASM_violin_plots_mark.pdf",w=10,h=4)
ggplot(df, aes(x=sample, y=ASM_score, colour = gender, fill=gender)) + facet_grid(~chrom) +
       geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(-1, 2)) +
       theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
       labs(x="Sample", y = "ASM Score")
dev.off()
