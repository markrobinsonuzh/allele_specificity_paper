## Figure 7 -- bsseq vs DAMEfinder
## Data: Adenoma data set
## rda files needed: "bsseq_dmrs_real_data.rda", "bsseq_real_data_samples.loci.rda", "real_odmrs_transformed.rda", "bsseq_samples_fit_real_data.rda"
## What we do: Plot the average methylation values per region for every region bsseq and DAMEfinder predict, and display them as normal vs adenoma.

library(ggplot2)
library(bsseq)

## Load data

# Load DMRs predicted by bsseq 
# (predicted with following coverage condition per position: at least 3 adenoma samples have coverage >=4 and at least 2 normal samples have coverage >=4)
load("../data/bsseq_dmrs_real_data.rda") # dmrs_subset

# Load BSseq object of the DMRs detected by bsseq 
load("../data/bsseq_real_data_samples.loci.rda") # samples.loci

# Load DAMEs predicted by DAMEfinder
load("../data/real_odmrs_transformed.rda") # odmrs

# Load BSseq object for the DAMEs. 
# This is the same bsseq object as for the DMRs of bsseq but before setting the condition mentioned above (so a bsseq object that has read in all the 'cov' files).
load("../data/bsseq_samples_fit_real_data.rda") # samples.fit

## getMeth

# The getMeth function return methylation estimates per region. 
# It estimates average methylation for a given region using information contained in a provided BSseq object.
# We use this function on the regions detected by bsseq and DAMEfinder on the adenoma data set.

region_meth_bsseq <- getMeth(samples.loci, regions = dmrs_subset, type = "smooth", what = "perRegion", confint = FALSE)
rownames(region_meth_bsseq) <- rownames(dmrs_subset)

region_meth_DAMEfinder <- getMeth(samples.fit, regions=odmrs, type = "smooth", what = "perRegion", confint = FALSE)
rownames(region_meth_DAMEfinder) <- rownames(odmrs)

## Get average methylation estimates for normal vs adenoma samples
x_avr_normal_meth <- c(rowMeans(region_meth_bsseq[,c(1,2,12)]), rowMeans(region_meth_DAMEfinder[,c(1,2,12)]))
y_avr_adenoma_meth <- c(rowMeans(region_meth_bsseq[,-c(1,2,12)]), rowMeans(region_meth_DAMEfinder[,-c(1,2,12)]))

## Get 'area' values for each detected region in bsseq and DAMEfinder
t_area <- c(abs(dmrs_subset$areaStat), odmrs$area)
method <- c(rep("bsseq", nrow(region_meth_bsseq)), rep("DAMEfinder", nrow(region_meth_DAMEfinder)))

plot_df <- data.frame(normal_meth=x_avr_normal_meth, adenoma_meth=y_avr_adenoma_meth, area=t_area, method=method)

## Plot All DAMEs and DMRs
# pdf("../figures/method_dmrs_meth_area_plot_2.pdf", w=10, h=8)
# qplot(normal_meth, adenoma_meth, data=plot_df, size=area, colour=factor(method), xlim = c(0,1), ylim = c(0,1), 
#       xlab = "Average Normal Methylation", ylab = "Average Adenoma Methylation")
# dev.off()

## Plot the top 50 DAMEs and top 50 DMRs
x_avr_normal_meth <- c(rowMeans(region_meth_bsseq[1:50,c(1,2,12)]), rowMeans(region_meth_DAMEfinder[1:50,c(1,2,12)]))
y_avr_adenoma_meth <- c(rowMeans(region_meth_bsseq[1:50,-c(1,2,12)]), rowMeans(region_meth_DAMEfinder[1:50,-c(1,2,12)]))

t_area <- c(abs(dmrs_subset[1:50,]$areaStat), odmrs[1:50,]$area)

method <- c(rep("bsseq", 50), rep("DAMEfinder", 50))

plot_df <- data.frame(normal_meth=x_avr_normal_meth, adenoma_meth=y_avr_adenoma_meth, area=t_area, method=method)

pdf("../figures/bsseq_DAMEfinder_comparison_top50.pdf", w=10, h=8)
qplot(normal_meth, adenoma_meth, data=plot_df, size=area, colour=factor(method), xlim = c(0,1), ylim = c(0,1), 
      xlab = "Average Normal Methylation", ylab = "Average Adenoma Methylation")
dev.off()

## Plot the top 1000 DAMEs and top 1000 DMRs
x_avr_normal_meth <- c(rowMeans(region_meth_bsseq[1:1000,c(1,2,12)]), rowMeans(region_meth_DAMEfinder[1:1000,c(1,2,12)]))
y_avr_adenoma_meth <- c(rowMeans(region_meth_bsseq[1:1000,-c(1,2,12)]), rowMeans(region_meth_DAMEfinder[1:1000,-c(1,2,12)]))

t_area <- c(abs(dmrs_subset[1:1000,]$areaStat), odmrs[1:1000,]$area)

method <- c(rep("bsseq", 1000), rep("DAMEfinder", 1000))

plot_df <- data.frame(normal_meth=x_avr_normal_meth, adenoma_meth=y_avr_adenoma_meth, area=t_area, method=method)

pdf("../figures/bsseq_DAMEfinder_comparison_top1000.pdf", w=10, h=8)
qplot(normal_meth, adenoma_meth, data=plot_df, size=area, colour=factor(method), xlim = c(0,1), ylim = c(0,1), 
      xlab = "Average Normal Methylation", ylab = "Average Adenoma Methylation")
dev.off()
