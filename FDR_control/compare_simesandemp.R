library(GenomicRanges)
library(DAMEfinder)
library(SummarizedExperiment)
library(ggplot2)

load("data/tupledames_cimp.RData")
sime_cimp <- GRanges(dames_cimp$chr, IRanges(dames_cimp$start, dames_cimp$end), 
                     pval = dames_cimp$pvalSimes,FDR = dames_cimp$FDR)
load("data/tupledames_noncimp.RData")
sime_non <- GRanges(dames_noncimp$chr, IRanges(dames_noncimp$start, dames_noncimp$end),
                    pval = dames_noncimp$pvalSimes,FDR = dames_noncimp$FDR)

#load("data/tupledames_cimp_emp.RData")
#load("data/tupledames_cimp_emp_repeat.RData")
load("data/tupledames_cimp_emp_repeat2.RData")
emp_cimp <- GRanges(dames_cimp$chr, IRanges(dames_cimp$start, dames_cimp$end),
                    pval = dames_cimp$pvalEmp,FDR = dames_cimp$FDR)
load("data/tupledames_noncimp_emp.RData")
emp_non <- GRanges(dames_noncimp$chr, IRanges(dames_noncimp$start, dames_noncimp$end),
                   pval = dames_noncimp$pvalEmp,FDR = dames_noncimp$FDR)

o <- findOverlaps(emp_cimp, sime_cimp)
tab <- data.frame(simes = c(-log10(sime_cimp$FDR[subjectHits(o)]),
                            -log10(sime_cimp$pval[subjectHits(o)])),
                  empirical = c(-log10(emp_cimp$FDR[queryHits(o)]),
                                -log10(emp_cimp$pval[queryHits(o)])),
                  val = factor(c(rep("FDR",length(o)), rep("pval",length(o))),
                               levels = c("pval", "FDR")))


a <- ggplot(tab) + 
  geom_bin2d(aes(simes, empirical)) + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  geom_abline(color = "darkgray") +
  facet_grid(~val) +
  theme_bw()

o <- findOverlaps(emp_non, sime_non)
tab <- data.frame(simes = c(-log10(sime_non$FDR[subjectHits(o)]),
                            -log10(sime_non$pval[subjectHits(o)])),
                  empirical = c(-log10(emp_non$FDR[queryHits(o)]),
                                -log10(emp_non$pval[queryHits(o)])),
                  val = factor(c(rep("FDR",length(o)), rep("pval",length(o))),
                               levels = c("pval", "FDR")))

b <- ggplot(tab) + 
  geom_bin2d(aes(simes, empirical)) + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  geom_abline(color = "darkgray") +
  facet_grid(~val) +
  theme_bw()
b
cowplot::plot_grid(a,b, nrow = 2, ncol = 1, labels = c("A", "B"))
ggsave("curvesNscatters/simesVsemp_both.png", width = 6, height = 6)


