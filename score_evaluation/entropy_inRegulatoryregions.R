
#Boxplots for Entropy, same as figure 3?
library(SummarizedExperiment)
library(ggplot2)
library(annotatr)
library(biomaRt)

load("data/malevsfem_ASM_mat.RData")

#Filter
ASM_mat <- ASM_mat[rowSums(
  !is.na(assays(ASM_mat)[["cov"]]) &
    assays(ASM_mat)[["cov"]] >= 10) == BiocGenerics::ncol(ASM_mat),] #3,198,678

#Add entropy
ME <- function(MM,UU,MU,UM,cov){
    if(MM != 0 & !is.na(MM)) pMM  <- (MM/cov) * log2(MM/cov)
    else if(is.na(MM)) pMM <- NA
    else pMM <- 0
    
    if(UU != 0 & !is.na(UU)) pUU  <- (UU/cov) * log2(UU/cov)
    else if(is.na(UU)) pUU <- NA
    else pUU <- 0
    
    if(MU != 0 & !is.na(MU)) pMU  <- (MU/cov) * log2(MU/cov)
    else if(is.na(MU)) pMU <- NA
    else pMU <- 0
    
    if(UM != 0 & !is.na(UM)) pUM  <- (UM/cov) * log2(UM/cov)
    else if(is.na(UM)) pUM <- NA
    else pUM <- 0
    
    return(0.5 * (-pMM - pUU - pMU - pUM))
}

ent <- sapply(1:6, function(number){
  mapply(ME,
         assay(ASM_mat,"MM")[,number],
         assay(ASM_mat,"UU")[,number],
         assay(ASM_mat,"MU")[,number],
         assay(ASM_mat,"UM")[,number],
         assay(ASM_mat,"cov")[,number])
})

assay(ASM_mat,"entropy") <- ent
save(ASM_mat, "data/malevsfem_ASM_mat_withent.RData")

#Prepare for plot
get_e <- function(mat){
  
  ASM_x <- mat[seqnames(mat) == "X",]
  ASM_3 <- mat[seqnames(mat) == "3",]
  
  asmX <- as.data.frame(assays(ASM_x)[["entropy"]])
  asm3 <- as.data.frame(assays(ASM_3)[["entropy"]])
  colnames(asmX) <- colnames(asm3) <- c(1:6)
  
  eX <- tidyr::gather(asmX, sample, valueasm)
  eX$chr <- "chrX"
  
  e3 <- tidyr::gather(asm3, sample, valueasm)
  e3$chr <- "chr3"
  
  e <- rbind(eX,e3)
  e$Gender <- 0
  e$Gender <- ifelse(e$sample %in% 1:3, "Male", "Female")
  return(e)
}

#### Replicate Figure 3 with entropies ####

## full
e <- get_e(ASM_mat)

myColor <- RColorBrewer::brewer.pal(9, "Set1")[3:4]
pfull <- ggplot(e) +
  geom_violin(aes(x=sample, y=valueasm, fill=Gender, color = Gender),
              trim = FALSE) +#, adjust = 1.5)+ #scale = "width") +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = "grey",
               alpha = 0,outlier.shape = NA) +
  #scale_y_continuous(trans='sqrt') +
  scale_fill_manual(values = myColor) +
  scale_color_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15),
        legend.position = "none") +
  labs(y="ASMtuple",x="") +
  facet_grid(~chr)


## proms
proms <- build_annotations(genome = 'hg38', annotations = "hg38_genes_promoters")
cpgsites <- rowRanges(ASM_mat) # 3,198,678
seqlevels(cpgsites) <- paste0("chr", seqlevels(cpgsites))

over <- findOverlaps(cpgsites, proms) #1,666,247
uniQ <- unique(queryHits(over))

ASM_mat.inproms <- ASM_mat[uniQ,]

e <- get_e(ASM_mat.inproms)
e$chr <- ifelse(e$chr == "chrX", "chrX (promoters)", e$chr)
e$chr <- ifelse(e$chr == "chr3", "chr3 (promoters)", e$chr)

pprom <- ggplot(e) +
  geom_violin(aes(x=sample, y=valueasm, fill=Gender, color = Gender),
              trim = FALSE ) +#, adjust = 1.5, scale = "width") +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = "grey",
               alpha = 0,outlier.shape = NA) +
  #scale_y_continuous(trans='sqrt') +
  scale_color_manual(values = myColor) +
  scale_fill_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15)) +
  labs(y="Entropy",x = "") +
  facet_grid(~chr)

## imprinted regs
impr <- read.csv("data/noplacenta.csv") #113

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
slopGR <- promoters(locsGR, upstream=1000, downstream=10)

cpgsites <- rowRanges(ASM_mat) 

#get ASM for imprinted genes
over <- findOverlaps(cpgsites, slopGR)
uniQ <- unique(queryHits(over))
ASM_mat.ingenes <- ASM_mat[uniQ,]

#get ASM for rest of genome
ASM_mat.rest <- ASM_mat[-uniQ,]

asmX <- as.data.frame(assays(ASM_mat.ingenes)[["entropy"]])
asm3 <- as.data.frame(assays(ASM_mat.rest)[["entropy"]])
colnames(asmX) <- colnames(asm3) <- c(1:6)

eX <- tidyr::gather(asmX, sample, valueasm)
eX$state <- "Imprinted"

e3 <- tidyr::gather(asm3, sample, valueasm)
e3$state <- "Genome"

e <- rbind(eX,e3)
e$Gender <- 0
e$Gender <- ifelse(e$sample %in% 1:3, "Male", "Female")

impr <- ggplot(e) +
  geom_violin(aes(x=sample, y=valueasm, fill=Gender, color = Gender),
              trim = FALSE) +#, adjust = 2) +# scale = "width") +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = "grey",
               alpha = 0,outlier.shape = NA) +
  #scale_y_continuous(trans='sqrt') +
  scale_fill_manual(values = myColor) +
  scale_color_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15),
        legend.position = "none") +
  ylab("Entropy") +
  facet_grid(~state)

p4 <- cowplot::plot_grid(pfull, pprom, impr, ncol=1, nrow = 3, 
                         labels = c("A","B","C"),
                         align = "v",
                         axis = "lr")
ggplot2::ggsave("curvesNscatters/chromboxplots_entropy.png", p4, 
                width = 10, height = 13)


#### Using annotation from Onuchic paper Figure1 ####
#I used annotar to get proms with-without cpgIslans, and phantom enhancers
#This is btw worst way of doing this

#in proms
proms <- build_annotations(genome = 'hg38', annotations = "hg38_genes_promoters")
isles <- build_annotations(genome = 'hg38', annotations = "hg38_cpg_islands")
over <- findOverlaps(proms, isles) 
proms_isles <- proms[unique(queryHits(over))]
#proms_non <- proms[-unique(queryHits(over))]

cpgsites <- rowRanges(ASM_mat) # 3,198,678
seqlevels(cpgsites) <- paste0("chr", seqlevels(cpgsites))
over <- findOverlaps(cpgsites, proms_isles) 
uniQ <- unique(queryHits(over))
ASM_mat.inproms_isles <- ASM_mat[uniQ,]
ASM_mat.inproms_non <- ASM_mat[-uniQ,]

#in isles
over <- findOverlaps(cpgsites, isles) 
uniQ <- unique(queryHits(over))
ASM_mat.inisles <- ASM_mat[uniQ,]

#in enhancers
enhan <- build_annotations(genome = 'hg38', annotations = "hg38_enhancers_fantom")
over <- findOverlaps(cpgsites, enhan)
uniQ <- unique(queryHits(over))
ASM_mat.inenhan <- ASM_mat[uniQ,]

#Join tables
asmP <- as.data.frame(assays(ASM_mat.inproms_isles)[["entropy"]])
asmPn <- as.data.frame(assays(ASM_mat.inproms_non)[["entropy"]])
asmI <- as.data.frame(assays(ASM_mat.inisles)[["entropy"]])
asmE <- as.data.frame(assays(ASM_mat.inenhan)[["entropy"]])
colnames(asmP) <- colnames(asmPn) <- colnames(asmI) <- colnames(asmE) <- c(1:6)

eP <- tidyr::gather(asmP, sample, valueasm)
eP$state <- "Promoters with CpG Islands"

ePn <- tidyr::gather(asmPn, sample, valueasm)
ePn$state <- "Promoters without CpG Islands"

eI <- tidyr::gather(asmI, sample, valueasm)
eI$state <- "CpG Islands"

eE <- tidyr::gather(asmE, sample, valueasm)
eE$state <- "Enhancers"


e <- rbind(eP,ePn,eI,eE)
e$Gender <- 0
e$Gender <- ifelse(e$sample %in% 1:3, "Male", "Female")

p4 <- ggplot(e) +
  geom_violin(aes(x=sample, y=valueasm, fill=Gender, color = Gender),
              trim = FALSE) +#, adjust = 2) +# scale = "width") +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = "grey",
               alpha = 0,outlier.shape = NA) +
  #scale_y_continuous(trans='sqrt') +
  scale_fill_manual(values = myColor) +
  scale_color_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15),
        legend.position = "none") +
  ylab("Entropy") +
  facet_wrap(~state)

ggplot2::ggsave("curvesNscatters/chromboxplots_entropy_inregsites.png", p4, 
                width = 10, height = 10)

#Same plot as above but for ASMtuple
asmP <- as.data.frame(assays(ASM_mat.inproms_isles)[["asm"]])
asmPn <- as.data.frame(assays(ASM_mat.inproms_non)[["asm"]])
asmI <- as.data.frame(assays(ASM_mat.inisles)[["asm"]])
asmE <- as.data.frame(assays(ASM_mat.inenhan)[["asm"]])
colnames(asmP) <- colnames(asmPn) <- colnames(asmI) <- colnames(asmE) <- c(1:6)

eP <- tidyr::gather(asmP, sample, valueasm)
eP$state <- "Promoters with CpG Islands"

ePn <- tidyr::gather(asmPn, sample, valueasm)
ePn$state <- "Promoters without CpG Islands"

eI <- tidyr::gather(asmI, sample, valueasm)
eI$state <- "CpG Islands"

eE <- tidyr::gather(asmE, sample, valueasm)
eE$state <- "Enhancers"


e <- rbind(eP,ePn,eI,eE)
e$ASMstatus <- ifelse(e$valueasm < 0.8, "noASM", "ASM")
e$Gender <- 0
e$Gender <- ifelse(e$sample %in% 1:3, "Male", "Female")

p5 <- ggplot(e) +
  geom_violin(aes(x=sample, y=valueasm, fill=Gender, color = Gender),
              trim = FALSE, adjust = 2) +# scale = "width") +
  geom_boxplot(aes(x=sample, y=valueasm, fill=Gender), color = "grey",
               alpha = 0,outlier.shape = NA) +
  scale_y_continuous(trans='sqrt') +
  scale_fill_manual(values = myColor) +
  scale_color_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15),
        legend.position = "none") +
  ylab("ASMtuple") +
  facet_wrap(~ASMstatus+state)
ggplot2::ggsave("curvesNscatters/chromboxplots_ASMtuple_inregsites.png", p5, 
                width = 10, height = 10)
ggplot2::ggsave("curvesNscatters/chromboxplots_ASMtuple_inregsites_byASMstate.png", p5, 
                width = 10, height = 10)