## Figure 7 -- bsseq vs DAMEfinder
## Data: Adenoma data set
## rda files needed: "bsseq_DAMEfinder_comparison_top50.rds", "bsseq_DAMEfinder_comparison_top1000.rds"
## extra file: "bsseq_DAMEfinder_comparison_all.rds" contains all the DAMEs and DMRs
## What we do: Plot the average methylation values per region for every region bsseq and DAMEfinder predict, and display them as normal vs adenoma.

## Plot top 50 DAMEs and top 50 DMRs
top50_df <- readRDS("../data/bsseq_DAMEfinder_comparison_top50.rds")

pdf("../figures/bsseq_DAMEfinder_comparison_top50.pdf", w=10, h=8)
qplot(normal_meth, adenoma_meth, data=top50_df, size=area, colour=factor(method), xlim = c(0,1), ylim = c(0,1), 
      xlab = "Average Normal Methylation", ylab = "Average Adenoma Methylation")
dev.off()

## Plot top 1000 DAMEs and top 1000 DMRs
top1000_df <- readRDS("../data/bsseq_DAMEfinder_comparison_top1000.rds")

pdf("../figures/bsseq_DAMEfinder_comparison_top1000.pdf", w=10, h=8)
qplot(normal_meth, adenoma_meth, data=top1000_df, size=area, colour=factor(method), xlim = c(0,1), ylim = c(0,1), 
      xlab = "Average Normal Methylation", ylab = "Average Adenoma Methylation")
dev.off()

