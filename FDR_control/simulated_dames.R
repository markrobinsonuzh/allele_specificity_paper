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
library(tidyr)

#### Set sim ####

load("data/derASM_fullCancer.RData")
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

#### play with clust length given maxGap####

#for sim1
# clust <- bumphunter::clusterMaker(as.character(seqnames(derASM)), start(derASM), maxGap = 20)
# max20 <- data.frame(clusL = rle(clust)$length, maxGap = 20)

#for sim2
clust <- bumphunter::clusterMaker(as.character(seqnames(derASM)), start(derASM), maxGap = 100)
maxcien <- data.frame(clusL = rle(clust)$length, maxGap = 100)

# clust <- bumphunter::clusterMaker(as.character(seqnames(derASM)), start(derASM), maxGap = 1000)
# maxmil <- data.frame(clusL = rle(clust)$length, maxGap = 1000)
# 
# clustab <- rbind(max20,maxcien,maxmil)
# ggplot(clustab, aes(clusL)) + geom_histogram() + facet_grid(~maxGap) + theme_bw()


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


#get real coordinates to start from
chr <- as.character(seqnames(derASM))
starts <- start(derASM)  
ends <- end(derASM) 

realregs <- data.frame(chr=sapply(cluster.ids,function(Index) chr[clust == Index][1]),
                       start=sapply(cluster.ids,function(Index) min(starts[clust == Index])),
                       end=sapply(cluster.ids, function(Index) max(ends[clust == Index])),
                       clusL=sapply(cluster.ids, function(Index) length(clust[clust == Index])))


#create 50 more simulations to run the methods

all_perf <- list()
all_points <- list()

for(j in 1:50){
  
  print(j)
  prop.clust <- x[,7:12]
  d <- qbeta(runif(length(diffClusts), minb, maxb), alpha, beta)
  #hist(d)


  ### Simulation ####

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

  #table(realregsGR$label)

  #head(prop.clust)
  #head(original)

  #re-do plots with added effects
  # means <- rowMeans(prop.clust)
  # var <- rowVars(prop.clust)
  # diffs <- apply(prop.clust, 1, function(w){mean(w[1:3]) - mean(w[4:6])})
  # dd <- as.data.frame(cbind(diffs, means,var))
  # ggplot(dd, aes(means, diffs)) + geom_point(alpha = 0.2) + theme_bw()
  # ggplot(dd, aes(means, var)) + geom_point(alpha = 0.2) + theme_bw()

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
  all_perf[[j]] <- cobraperf@fdrtprcurve
  all_points[[j]] <- cobraperf@fdrtpr
}

#### set up to plot all sims ####
#lines
tpr <- lapply(all_perf, function(x){x$TPR})

allperftab <- data.frame(sim = rep(1:50, lengths(tpr)),
                         FDR = unlist(lapply(all_perf, function(x){x$FDR})),
                         TPR = unlist(tpr),
                         method = unlist(lapply(all_perf, function(x){x$method})))

allperftab <- unite(allperftab, unique_id, c(sim, method), sep="_", remove = FALSE)

#allperftabfilt <- allperftab[!is.nan(allperftab$FDR) & !is.infinite(allperftab$FDR),]

# summperf <- allperftab %>%
#   filter(!is.na(FDR)) %>%
#   dplyr::group_by(method, sim) %>%
#   dplyr::summarise(meanTPR=mean(TPR), meanFDR = mean(FDR)) %>%
#   as.data.frame()

#points
tpr <- lapply(all_points, function(x){x$TPR})

allpointtab <- data.frame(sim = rep(1:50, lengths(tpr)),
                          FDR = unlist(lapply(all_points, function(x){x$FDR})),
                          TPR = unlist(tpr),
                          method = unlist(lapply(all_points, function(x){x$method})),
                          thr = unlist(lapply(all_points, function(x){x$thr})),
                          satis = unlist(lapply(all_points, function(x){x$satis})))

summpoints <- allpointtab %>%
  dplyr::group_by(method, thr) %>%
  dplyr::summarise(meanTPR=mean(TPR), meanFDR = mean(FDR)) %>%
  as.data.frame()

summpoints$thr <- as.numeric(gsub("thr","",summpoints$thr))
summpoints$satis <- ifelse(summpoints$meanFDR <= summpoints$thr,16,21)

myColor <- RColorBrewer::brewer.pal(8, "Set1")

ggplot(allperftab) +
  geom_line(aes(FDR, TPR, color=method, group=unique_id), alpha = 0.11) +
  #geom_smooth(aes(FDR, TPR, color=method, group=unique_id), na.rm = TRUE, method = "auto") +
  #geom_line(data = summperf, aes(meanFDR, meanTPR, color=method, group=unique_id)) +
  scale_x_continuous(trans='sqrt') +
  scale_color_manual(values = myColor) +
  labs(color = "Method") +
  geom_vline(xintercept = c(0.01,0.05,0.1), linetype = 2) +
  geom_line(data = summpoints, aes(x = meanFDR, y = meanTPR,color=method), size = 1) +
  geom_point(data = summpoints, aes(x = meanFDR,y = meanTPR,color=method, shape = satis), 
             size = 5, fill = "white") +
  scale_shape_identity() +
  theme_bw()
ggsave("curvesNscatters/powerFDR_means.png")







