
#tuple cimp simes
load("data/tupledames_cimp.RData")
dames_cimp_tuple <- dames_cimp[dames_cimp$FDR<0.05,] #4051
dames_cimp_tuple$state <- ifelse(dames_cimp_tuple$meanTstat > 0, "up","down")
tupl <- GRanges(dames_cimp_tuple$chr, IRanges(dames_cimp_tuple$start, dames_cimp_tuple$end), 
               pval = dames_cimp_tuple$FDR, )

#snp cimp simes
dames_cimp_snp <- find_dames(derASM, mod, coef = 1, contrast = cont, maxGap = 100) #2219
dames_cimp_snp
dames_cimp_snp$state <- ifelse(dames_cimp_snp$meanTstat > 0, "up","down")
snp <- GRanges(dames_cimp_snp$chr, IRanges(dames_cimp_snp$start, dames_cimp_snp$end), 
                pval = dames_cimp_snp$FDR)

over <- findOverlaps(tupl, snp, minoverlap = 1) #47
tupl_sub <- dames_cimp_tuple[queryHits(over),]
snp_sub <- dames_cimp_snp[subjectHits(over),]

x <- which(snp_sub$state != tupl_sub$state)

tupl_sub[x,]
snp_sub[x,]
