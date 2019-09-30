#!/usr/bin/env Rscript
# chmod +x
# run as [R < scriptName.R --no-save]

#########################################################################################
# Plot Null and real areas from empirical correction, from single simulation
#
# TBS-seq
#
#
# Stephany Orjuela, September 2019
#########################################################################################


empirical_pval2 <- function(presa, design, coeff, smooth,
                           maxPerms = 10, Q, maxGap, method, ...){

  sampleSize <- table(design[, coeff])

  #Since we allow the user to choose the coeff from the design.matrix, I will
  #never have more than 2 coefs, as dmrseq does

  if(length(unique(design[, coeff])) == 2 &&

     #limit possible number of perms!
     choose(nrow(design), min(sampleSize)) < 5e5 ){

    perms <- utils::combn(seq(1, nrow(design)), min(sampleSize))

    # Remove redundant permutations (if balanced)
    if(length(unique(table(design[,coeff]))) == 1){
      perms <- perms[, seq_len(ncol(perms)/2)]
    }

    # restrict to unique permutations that don't include any
    # groups consisting of all identical conditions
    rmv <- NULL
    for(p in seq_len(ncol(perms))){
      if(length(unique(design[perms[,p],coeff])) == 1){
        rmv <- c(rmv, p)
      }
    }

    if(length(rmv) > 0 ) perms <- perms[,-rmv]

    # subsample permutations based on similarity to original partition
    # gives preference to those with the least similarity
    if(maxPerms < ncol(perms)){
      similarity <- apply(perms, 2, function(x) {
        max(table(design[x,coeff]))
      })
      perms.all <- perms
      perms <- NULL
      levs <- sort(unique(similarity))
      l <- 1
      num <- 0
      while(!(num == maxPerms) && l <= length(levs)) {
        keep <- sample(which(similarity == levs[l]),
                       min(maxPerms-num, sum(similarity == levs[l])) )
        perms <- cbind(perms, perms.all[,keep])
        l <- l + 1
        num <- ncol(perms)
      }
    }
  } else message("Too many samples!")

  #Detect permuted dames
  message("Generating ", ncol(perms), " permutations", appendLF = TRUE)
  areas <- apply(perms, 2, function(i){

    reorder <- i
    designr <- design

    if(length(unique(design[, coeff])) == 2 &&
       !nrow(perms) == nrow(designr)){
      designr[,coeff] <- 0
      designr[reorder, coeff] <- 1
    } else {
      designr[, coeff] <- designr[reorder, coeff]
    }

    sa_perm <- get_tstats(presa,
                          design = designr,
                          coef = coeff,
                          maxGap = maxGap,
                          smooth  = smooth,
                          verbose = FALSE,
                          method = method,
                          ...)

    # choose smoothed if true
    if(smooth){
      sm_tstat <- S4Vectors::mcols(sa_perm)$smooth_tstat
    } else {
      sm_tstat <- S4Vectors::mcols(sa_perm)$tstat
    }

    #choose position to find regions
    if(names(assays(sa_perm))[1] == "asm"){
      midpt <- S4Vectors::mcols(sa_perm)$midpt
    } else {
      midpt <- BiocGenerics::start(sa_perm)
    }

    K <- stats::quantile(abs(sm_tstat), Q, na.rm=TRUE)
    permrf <- bumphunter::regionFinder(
      x = sm_tstat,
      chr = as.character(GenomeInfoDb::seqnames(sa_perm)),
      pos = midpt,
      cluster = S4Vectors::mcols(sa_perm)$cluster,
      cutoff = K, #if use same K always, doesnt generate any regions
      maxGap = maxGap,
      verbose = FALSE)

    return(abs(permrf$area))

  })

  #if( put a condition to check if areas actually has something
  all_areas <- sort(unlist(areas))
  total_areas <- length(all_areas)
  return(all_areas)
}


load("data/derASM_fullCancer.RData")
derASM <- GenomeInfoDb::sortSeqlevels(derASM) #only necessary for old calc_derivedasm()
derASM <- sort(derASM)

#use only the sites completely covered by all samples
filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) >= 12
derASM <- derASM[filt,]
x <- assay(derASM,"der.ASM")

##Get only norm samples
prop.clust <- x[,7:12]
original <- prop.clust

clust <- bumphunter::clusterMaker(as.character(seqnames(derASM)), start(derASM), maxGap = 100)


set.seed(20) # params for very obvious regions
alpha <- 1
beta <- 2.5
minb <- 0.35 # 0.15 too small for lmfit to consider it a difference
maxb <- 0.75 

pDiff <- 0.2 #this should affect the k choice
cluster.ids <- unique(clust)
diffClusts <- 1:floor(pDiff*length(cluster.ids))


#get real coordinates to start from
chr <- as.character(seqnames(derASM))
starts <- start(derASM)  
ends <- end(derASM) 

realregs <- data.frame(chr=sapply(cluster.ids,function(Index) chr[clust == Index][1]),
                       start=sapply(cluster.ids,function(Index) min(starts[clust == Index])),
                       end=sapply(cluster.ids, function(Index) max(ends[clust == Index])),
                       clusL=sapply(cluster.ids, function(Index) length(clust[clust == Index])))

prop.clust <- x[,7:12]
d <- qbeta(runif(length(diffClusts), minb, maxb), alpha, beta)



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
realregsGR <- realregsGR[filt] 

#build a sumExp with new data
fakeDerAsm <- derASM[,7:12]
assay(fakeDerAsm, "der.ASM") <- prop.clust
grp <- factor(c(rep("CRC",3),rep("NORM",3)), levels = c("NORM", "CRC"))
mod <- model.matrix(~grp)


#Get null areas
x <- empirical_pval2(presa = fakeDerAsm, design = mod, coeff = 2, smooth = TRUE,
                     maxPerms = 10, Q = 0.2, maxGap = 100, method = "ls")



#Get real areas
sat <- get_tstats(sa = fakeDerAsm, design = mod,
                  maxGap = 100,
                  coef = 2
                  )

sm_tstat <- S4Vectors::mcols(sat)$tstat

#choose position to find regions
midpt <- BiocGenerics::start(sat)

#detect dames
K <- stats::quantile(abs(sm_tstat), 0.2, na.rm=TRUE)

regs2 <- bumphunter::regionFinder(
    x = sm_tstat,
    chr = as.character(GenomeInfoDb::seqnames(sat)),
    pos = midpt,
    cluster = S4Vectors::mcols(sat)$cluster,
    cutoff = K,
    maxGap = 100,
    assumeSorted = TRUE,
    order = FALSE,
    verbose = TRUE)



#Set table to plot

nulls <- data.frame(areas = x, type = "nulls")
real <- data.frame(areas = regs2$area, type = "obtained")
toplot <- rbind(nulls, real)

ggplot(nulls, aes(areas)) + geom_histogram()


ggplot(toplot) +
  geom_histogram(aes(areas, fill = type, color = type), bins = 50, alpha = 0.2) +
  scale_x_continuous(trans = "log2") +
  theme_bw()

ggsave("curvesNscatters/areas_nullandreal_02_1sim.png")

#extra check
all_areas <- sort(unlist(x))
total_areas <- length(x)

pvalEmp <- sapply(regs2$area, function(a){
  pperm <- (sum(all_areas > abs(a)) + 1) / (total_areas + 1)
  return(pperm)
})
