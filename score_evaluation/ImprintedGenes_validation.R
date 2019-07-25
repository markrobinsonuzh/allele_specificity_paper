#########################################################################################
# Get ASM scores for known imprinted regions, and compare to rest of genome
#
# WGBS-seq data females Vs Males, BluePrint project
# NOTE: Files included here too big to include in repo
#
# Stephany Orjuela, July 2019
#########################################################################################

library(biomaRt)
library(GenomicRanges)
library(SummarizedExperiment)
library(ggplot2)


#load known imprinted regions
impr <- read.csv("noplacenta.csv") #113

#get coords
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
imprlocs <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position',
                               'strand', "hgnc_symbol", "external_gene_name",
                               "entrezgene_trans_name"),
      filters='hgnc_symbol',
      values=impr$Gene,
      mart=ensembl) #89

locsGR <- GRanges(imprlocs$chromosome_name, IRanges(imprlocs$start_position,
                                                    imprlocs$end_position),
                  strand = ifelse(imprlocs$strand == 1, "+", "-"),
                  hgnc_symbol = imprlocs$hgnc_symbol) 

#get promoters from these genes
slopGR <- promoters(locsGR, upstream=2000, downstream=10)

#or promoters from annotar (all tx)
#ind <- which(proms$symbol %in% impr$Gene)
#slopGR <- proms[ind,] #1330, 2035
 

#get full ASM mat
load("malevsfem_ASM_mat.RData")

#filter

cpgsites <- rowRanges(ASM_mat) # 3,198,678
#seqlevels(cpgsites) <- paste0("chr", seqlevels(cpgsites))

#get ASM for imprinted genes
over <- findOverlaps(cpgsites, slopGR)
uniQ <- unique(queryHits(over))
ASM_mat.ingenes <- ASM_mat[uniQ,]

#get ASM for rest of genome
ASM_mat.rest <- ASM_mat[-uniQ,]

asmX <- as.data.frame(assays(ASM_mat.ingenes)[["asm"]])
asm3 <- as.data.frame(assays(ASM_mat.rest)[["asm"]])
colnames(asmX) <- colnames(asm3) <- c(1:6)

eX <- tidyr::gather(asmX, sample, valueasm)
eX$state <- "Imprinted"
#eX$Gender <- 0
#eX$Gender <- ifelse(eX$sample %in% 1:3, "Male", "Female")


e3 <- tidyr::gather(asm3, sample, valueasm)
e3$state <- "Genome"

e <- rbind(eX,e3)
e$Gender <- 0
e$Gender <- ifelse(e$sample %in% 1:3, "Male", "Female")

#plot
myColor <- RColorBrewer::brewer.pal(9, "Set1")[3:4]
greys <- RColorBrewer::brewer.pal(9, "Greys")[7]

impr <- ggplot(e) +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = greys) +
  scale_fill_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15)) +
  ylab("ASMtuple") +
  facet_grid(~state)

ggsave("curvesNscatters/imprintboxplot_promgene.png",impr, width = 10, height = 10)

