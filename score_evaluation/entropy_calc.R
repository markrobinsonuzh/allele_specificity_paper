#Calculate shannon entropy from tuple methylation

load("data/tupleASM_fullCancer.RData") 


library(GenomicRanges)
library(SummarizedExperiment)
library(ggplot2)

number <- 1

#full entropy, from Guo 2017

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


#Calculate entropy
H1 <- mapply(ME,
             assay(ASM,"MM")[,number],
             assay(ASM,"UU")[,number],
             assay(ASM,"MU")[,number],
             assay(ASM,"UM")[,number],
             assay(ASM,"cov")[,number])

#ASM
asm <- assay(ASM,"asm")

ggplot() +
  geom_bin2d(aes(H1,asm[,number])) + 
  scale_fill_distiller(palette='RdBu', trans='log10') +
  theme_bw()

## Compare entropy between CpGs with SD-ASM and not with SD-ASM
load("data/derASM_fullcancer2.RData") #single site derived-true asm
load("data/tuplederASM_fullCancer2.RData")

#Get CpGs that have info of link to a het SNP
keyGR <- rowRanges(derASM) #1,453,873
start(keyGR) <- end(keyGR) <- start(keyGR) - 1
rm(derASM)

asmscoreGR <- rowRanges(ASM) #3,589,472
over <- findOverlaps(asmscoreGR, keyGR)
ASM <- ASM[unique(sort(queryHits(over))),] #957,310

asmscoreGR <- rowRanges(ASM)

#Add scores to GR
asmscoreGR$ASMsnp <- abs(assays(tuple.derASM)[["der.ASM"]][,number])

asmscoreGR$ASMtuple <- 0
asmscoreGR$ASMtuple <- assays(ASM)[["asm"]][,number]
asmscoreGR$ASMtuple <- ifelse(is.na(asmscoreGR$ASMtuple),0,asmscoreGR$ASMtuple)

#Get H1

asmscoreGR$Entropy <- 0
H1 <- mapply(ME,
             assay(ASM,"MM")[,number],
             assay(ASM,"UU")[,number],
             assay(ASM,"MU")[,number],
             assay(ASM,"UM")[,number],
             assay(ASM,"cov")[,number])
asmscoreGR$Entropy <- H1

#Remove things without SNPs, cov below 10
asmscoreGR <- asmscoreGR[!is.na(assay(tuple.derASM,"snp.table")[,number]) &
                         !is.na(assay(tuple.derASM,"ref.cov")[,number] + assay(tuple.derASM,"alt.cov")[,number]) &
                         (assay(tuple.derASM,"ref.cov")[,number] + assay(tuple.derASM,"alt.cov")[,number]) >= 10 
                         ] #only SNP:279,191, all:265,408

#plot
x <- mcols(asmscoreGR)
x$State <- ifelse(x$ASMsnp >= 0.8, "SD-ASM", "no SD-ASM")
x <- as.data.frame(x)

table(x$State)

p1 <- ggplot(x) +
  geom_bin2d(aes(Entropy, ASMtuple)) + 
  scale_fill_distiller(palette='RdBu', trans='log10') +
  labs(x = "ME methylation entropy")+
  theme_bw()

ggplot2::ggsave(filename = "curvesNscatters/EntropyVsASMtuple_inSD.png",p1)


ggplot(x) +
  geom_density(aes(ASMsnp, color = State, fill = State), alpha = 0.2) +
    theme_bw() 

p1 <- ggplot(x) +
  geom_density(aes(Entropy, color = State, fill = State), alpha = 0.2) +
  theme_bw()

p2 <- ggplot(x) +
  geom_density(aes(ASMtuple, color = State, fill = State), alpha = 0.2) +
  #geom_histogram(aes(ASMtuple, color = State, fill = State), alpha = 0.2) +
  scale_y_continuous(trans='sqrt') +
  theme_bw()

p3 <- cowplot::plot_grid(p1,p2, nrow = 1, ncol = 2,
                         labels = c("A","B"))

ggplot2::ggsave(filename = "curvesNscatters/EntropyandASMtuple_dists.png",p3, width = 10, height = 4)


