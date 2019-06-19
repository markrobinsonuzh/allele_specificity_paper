#!/usr/bin/env Rscript
# chmod +x
# run as [R < scriptName.R --no-save]

#########################################################################################
# Benchmark p-val assigment strategies with simulations 
#
# TBS-seq data CRCs Vs Norm
#
#
# Stephany Orjuela, May 2019
#########################################################################################

library(SummarizedExperiment)
library(ggplot2)
library(iCOBRA)

#### Set sim ####

load("../derASM_fullCancer.RData")
derASM <- GenomeInfoDb::sortSeqlevels(derASM) #only necessary for old calc_derivedasm()
derASM <- sort(derASM)

#use only the sites completely covered by all samples
filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) >= 12
derASM <- derASM[filt,] #9073
x <- assay(derASM,"der.ASM")

##Get only norm samples
prop.clust <- x[,7:12]
original <- prop.clust

means <- rowMeans(prop.clust)
diffs <- apply(prop.clust, 1, function(w){mean(w[1:3]) - mean(w[4:6])})
var <- rowVars(prop.clust)
dd <- as.data.frame(cbind(var, means, diffs))
head(dd)

#MD plot
ggplot(dd, aes(means, diffs)) + geom_point(alpha = 0.2) +
  theme_bw()

#MV plot
ggplot(dd, aes(means, var)) + geom_point(alpha = 0.2) + theme_bw()

# set true regions
#### play with clust length given maxGap####

#for sim1
clust <- bumphunter::clusterMaker(as.character(seqnames(derASM)), start(derASM), maxGap = 20)
max20 <- data.frame(clusL = rle(clust)$length, maxGap = 20)

#for sim2
clust <- bumphunter::clusterMaker(as.character(seqnames(derASM)), start(derASM), maxGap = 100)
maxcien <- data.frame(clusL = rle(clust)$length, maxGap = 100)

clust <- bumphunter::clusterMaker(as.character(seqnames(derASM)), start(derASM), maxGap = 1000)
maxmil <- data.frame(clusL = rle(clust)$length, maxGap = 1000)

clustab <- rbind(max20,maxcien,maxmil)
ggplot(clustab, aes(clusL)) + geom_histogram() + facet_grid(~maxGap) + theme_bw()


#### inverse sampling ####
#inverse sampling with truncated beta
set.seed(20) # params for very obvious regions
alpha <- 1
beta <- 2.5
minb <- 0.35 # 0.15 too small for lmfit to consider it a difference
maxb <- 0.75 

pDiff <- 0.5 #this should affect the k choice
cluster.ids <- unique(clust) #3229, 1038
diffClusts <- 1:floor(pDiff*length(cluster.ids)) #645
d <- qbeta(runif(length(diffClusts), minb, maxb), alpha, beta)
hist(d)

### Simulation 1 ####

chr <- as.character(seqnames(derASM))
starts <- start(derASM)
ends <- end(derASM)

realregs <- data.frame(chr=sapply(cluster.ids,function(Index) chr[clust == Index][1]),
                       start=sapply(cluster.ids,function(Index) min(starts[clust == Index])),
                       end=sapply(cluster.ids, function(Index) max(ends[clust == Index])),
                       clusL=sapply(cluster.ids, function(Index) length(clust[clust == Index])))

#make real GRanges
realregsGR <- GRanges(realregs$chr, IRanges(realregs$start, realregs$end), 
                      clusL = realregs$clusL,
                      label = c(rep(1,length(diffClusts)), 
                                rep(0,(length(cluster.ids)-length(diffClusts)))))

filt <- realregsGR$clusL != 1
realregsGR <- realregsGR[filt]

table(realregsGR$label) #1767

#Go through the clusters that should be DAMES and add effect size
for(i in diffClusts){
  
  #print(i)
  #get CpGs per cluster
  cpgs <- which(clust == cluster.ids[i])
  
  #print(cpgs)
  #randomly choose which group is diff
  ran <- sample(c(1,2),1)
  if(ran == 1) {group <- 1:3} else {group <- 4:6}
  
  #get cluster ASMsnp mean (if more than one sample)
  if(length(cpgs) > 1){
    DMRmean <- mean(rowMeans(prop.clust[cpgs,]))
  } else{
    DMRmean <- mean(prop.clust[cpgs,])
  }
  
  #sign is deterministic: 
  #if the DMR mean (across samples and loci) is below
  #effect size 0.5, sign is positive
  
  if(DMRmean < 0.5) {sign <- 1} else {sign <- -1}
  
  #if any of the values goes outside of [0,1], keep the original prop (second)
  prop.clust[cpgs,group] <- original[cpgs,group] + (d[i] * sign)
  
  if(any(prop.clust[cpgs,group] < 0 | prop.clust[cpgs,group] > 1)){ 
    w <- which(prop.clust[cpgs,group] < 0 | prop.clust[cpgs,group] > 1)
    prop.clust[cpgs,group][w] <- original[cpgs,group][w]
  }
  
}

head(prop.clust)
head(original)

#re-do plots with added effects
means <- rowMeans(prop.clust)
var <- rowVars(prop.clust)
diffs <- apply(prop.clust, 1, function(w){mean(w[1:3]) - mean(w[4:6])})
dd <- as.data.frame(cbind(diffs, means,var))
ggplot(dd, aes(means, diffs)) + geom_point(alpha = 0.2) + theme_bw()
ggplot(dd, aes(means, var)) + geom_point(alpha = 0.2) + theme_bw()

#build a sumExp with new data
fakeDerAsm <- derASM[,7:12]
assay(fakeDerAsm, "der.ASM") <- prop.clust
grp <- factor(c(rep("CRC",3),rep("NORM",3)), levels = c("NORM", "CRC"))
mod <- model.matrix(~grp)
#tstat <- get_tstats(fakeDerAsm, mod, maxGap = 20)

### Simulation 2 ####

chr <- as.character(seqnames(derASM))
starts <- start(derASM)
ends <- end(derASM)

realregs <- data.frame(chr=sapply(cluster.ids,function(Index) chr[clust == Index][1]),
                       start=sapply(cluster.ids,function(Index) min(starts[clust == Index])),
                       end=sapply(cluster.ids, function(Index) max(ends[clust == Index])),
                       clusL=sapply(cluster.ids, function(Index) length(clust[clust == Index])))


for(i in diffClusts){
  
  #get CpGs per cluster that will have spike-in
  cpgs <- which(clust == cluster.ids[i])
  
  #choose number of CpGs diff per regions, and from what position
  if(length(cpgs) > 1){
    numdiff <- sample(1:length(cpgs), 1)
    maxpos <- length(cpgs) - numdiff + 1
    posdiff <- sample(1:maxpos,1)
    cpgs <- cpgs[1:posdiff]
    
    #reset region start end ends
    realregs$start[i] <- min(starts[cpgs])
    realregs$end[i] <- max(ends[cpgs])
    realregs$clusL[i] <- length(cpgs)
  }
  
  
  
  #randomly choose which group is diff
  ran <- sample(c(1,2),1)
  if(ran == 1) {group <- 1:3} else {group <- 4:6}
  
  #get cluster ASMsnp mean (if more than one sample)
  if(length(cpgs) > 1){
    DMRmean <- mean(rowMeans(prop.clust[cpgs,]))
  } else{
    DMRmean <- mean(prop.clust[cpgs,])
  }
  
  #sign is deterministic: 
  #if the DMR mean (across samples and loci) is below
  #effect size 0.5, sign is positive
  
  if(DMRmean < 0.5) {sign <- 1} else {sign <- -1}
  
  #if any of the values goes outside of [0,1], keep the original prop (second)
  prop.clust[cpgs,group] <- original[cpgs,group] + (d[i] * sign)
  
  if(any(prop.clust[cpgs,group] < 0 | prop.clust[cpgs,group] > 1)){ 
    w <- which(prop.clust[cpgs,group] < 0 | prop.clust[cpgs,group] > 1)
    prop.clust[cpgs,group][w] <- original[cpgs,group][w]
  }
  
}

#make real GRanges
realregsGR <- GRanges(realregs$chr, IRanges(realregs$start, realregs$end), 
                      clusL = realregs$clusL,
                      label = c(rep(1,length(diffClusts)), 
                                rep(0,(length(cluster.ids)-length(diffClusts)))))


filt <- realregsGR$clusL != 1
realregsGR <- realregsGR[filt] #773

table(realregsGR$label)

head(prop.clust)
head(original)

#re-do plots with added effects
means <- rowMeans(prop.clust)
var <- rowVars(prop.clust)
diffs <- apply(prop.clust, 1, function(w){mean(w[1:3]) - mean(w[4:6])})
dd <- as.data.frame(cbind(diffs, means,var))
ggplot(dd, aes(means, diffs)) + geom_point(alpha = 0.2) + theme_bw()
ggplot(dd, aes(means, var)) + geom_point(alpha = 0.2) + theme_bw()

#build a sumExp with new data
fakeDerAsm <- derASM[,7:12]
assay(fakeDerAsm, "der.ASM") <- prop.clust
grp <- factor(c(rep("CRC",3),rep("NORM",3)), levels = c("NORM", "CRC"))
mod <- model.matrix(~grp)

#### Apply all methods ####

#simes
regs <- find_dames(fakeDerAsm, mod, maxGap = 100)
regsGR <- GRanges(regs$chr, IRanges(regs$start, regs$end), 
                  clusterL = regs$clusterL, pval = regs$pvalSimes, FDR = regs$FDR)

# regs <- find_dames(fakeDerAsm, mod, maxGap = 20, method = "robust") #3229
# regsRobGR <- GRanges(regs$chr, IRanges(regs$start, regs$end), 
#                   clusterL = regs$clusterL, pval = regs$pvalSimes, FDR = regs$FDR)

#test trend
# regs <- find_dames(fakeDerAsm, mod, maxGap = 20, trend = TRUE)
# regsRobTrendGR <- GRanges(regs$chr, IRanges(regs$start, regs$end), 
#                      clusterL = regs$clusterL, pval = regs$pvalSimes, FDR = regs$FDR)


#empirical

regs2 <- find_dames(fakeDerAsm, mod, maxGap = 100, pvalAssign = "empirical", Q = 0.2)
regs1GR <- GRanges(regs2$chr, IRanges(regs2$start, regs2$end), segmentL = regs2$segmentL, 
                   clusterL = regs2$clusterL, pval  = regs2$pvalEmp, FDR = regs2$FDR)

regs2 <- find_dames(fakeDerAsm, mod, maxGap = 100, pvalAssign = "empirical", Q = 0.5)
regs2GR <- GRanges(regs2$chr, IRanges(regs2$start, regs2$end), segmentL = regs2$segmentL, 
                  clusterL = regs2$clusterL, pval  = regs2$pvalEmp, FDR = regs2$FDR)

regs2 <- find_dames(fakeDerAsm, mod, maxGap = 100, pvalAssign = "empirical", Q = 0.8)
regs3GR <- GRanges(regs2$chr, IRanges(regs2$start, regs2$end), segmentL = regs2$segmentL, 
                   clusterL = regs2$clusterL, pval  = regs2$pvalEmp, FDR = regs2$FDR)


#### build tables with pval methods ####
# pvalmat <- data.frame(matrix(1, nrow = length(cluster.ids), ncol = 5))
# fdrmat <- data.frame(matrix(1, nrow = length(cluster.ids), ncol = 5))
# colnames(pvalmat) <- colnames(fdrmat) <- c("simes","simes.robust",# "simes.trend",
#                                            "perms","perms.robust","bumpperms")

pvalmat <- data.frame(matrix(1, nrow = length(realregsGR), ncol = 4))
fdrmat <- data.frame(matrix(1, nrow = length(realregsGR), ncol = 4))
colnames(pvalmat) <- colnames(fdrmat) <- c("simes",
                                           "perms_02",
                                           "perms_05",
                                           "perms_08")

#simes
over <- findOverlaps(realregsGR, regsGR, type = "within")
pvalmat$simes[queryHits(over)] <- mcols(regsGR)$pval[subjectHits(over)]
fdrmat$simes[queryHits(over)] <- mcols(regsGR)$FDR[subjectHits(over)]

#perms.0.2
over <- findOverlaps(realregsGR, regs1GR, type = "within")
pvalmat$perms_02[queryHits(over)] <- mcols(regs1GR)$pval[subjectHits(over)]
fdrmat$perms_02[queryHits(over)] <- mcols(regs1GR)$FDR[subjectHits(over)]

#perms.0.5
over <- findOverlaps(realregsGR, regs2GR, type = "within")
pvalmat$perms_05[queryHits(over)] <- mcols(regs2GR)$pval[subjectHits(over)]
fdrmat$perms_05[queryHits(over)] <- mcols(regs2GR)$FDR[subjectHits(over)]

#perm.0.8
over <- findOverlaps(realregsGR, regs3GR, type = "within")
pvalmat$perms_08[queryHits(over)] <- mcols(regs3GR)$pval[subjectHits(over)]
fdrmat$perms_08[queryHits(over)] <- mcols(regs3GR)$FDR[subjectHits(over)]


#### plot powerFDR ####
#generate truth + facet table
truth <- as.data.frame(mcols(realregsGR))
#change clusL to num.CpGs

#run iCOBRa
cobradat <- COBRAData(pval = pvalmat,
                      padj = fdrmat,
                      truth = truth)

#single plot
cobraperf <- calculate_performance(cobradat, binary_truth = "label", 
                                   cont_truth = "label",
                                   aspects = c("fdrtpr","fdrtprcurve"),
                                   thrs = c(0.01, 0.05, 0.1))

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Set1",facetted = FALSE)

#title = "diff.regs = 0.5, sim = 1, maxGap = 20")
p1 <- plot_fdrtprcurve(cobraplot) +
  scale_x_continuous(trans='sqrt') +
  theme(axis.text.x = element_text(size = 8,angle=0),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        legend.position="none")  


##accuracy vs numCpGs (clusL)
cobraperf <- calculate_performance(cobradat, binary_truth = "label", 
                                   cont_truth = "label", splv = "clusL",
                                   aspects = c("fdrtpr","fdrtprcurve"), 
                                   maxsplit = Inf,
                                   #thrs = c(0.01, 0.05, 0.1))
                                   thrs = 0.05)

x <- cobraperf@fdrtpr
x$accuracy <- (x$TP+x$TN)/(x$TP+x$TN+x$FP+x$FN)
x$numCGs <- as.numeric(gsub("(perms_0[2-8]|simes)_(clusL:|overall)","",x$fullmethod))
x <- tidyr::drop_na(x)
myColor <- RColorBrewer::brewer.pal(9, "Set1")

p2 <- ggplot(x, aes(x = numCGs, y = accuracy, color = method)) + #shape = thr)) + 
  geom_line() +
  geom_point() +
  scale_color_manual(values = myColor) +
  scale_shape_manual(values = 1:3) +
  theme_bw()

cowplot::plot_grid(p1, p2, ncol=2, nrow = 1, align="h", rel_widths = c(1,1.5), 
                         rel_heights = 1)

