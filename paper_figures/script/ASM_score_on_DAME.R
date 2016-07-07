#  Here we plot the ASM score and the smoothed values on a specific DAME on chromosome 8: 141108266-141111068
# Figure 10

# load asm_matrix
load("../data/asm_matrix.rda")

library(bumphunter)

dame_chr <- "chr8"
dame_start <- 141108266 
dame_end <- 141111068

# the rownames of the matrix are called as follows: chr.pos1.pos2
rows <- rownames(sm)
asm_chr <- limma::strsplit2(rows, ".", fixed=T)[,1]
asm_pos1 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,2])
asm_pos2 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,3])

# take the midpoint of pos1 and pos2 as the genomic position for each row
midpt <- floor((asm_pos2 - asm_pos1)/2)
asm_pos <- asm_pos1 + midpt

# mke sure the positions are ordered 
o <- order(asm_pos)
asm_matrix <- asm_matrix[o,]
asm_pos <- asm_pos[o]
asm_chr <- asm_chr[o]

# the column names in the matrix indicate the type and number of each sample
asm_cols <- colnames(asm_matrix)

# get smoothed values of the ASM scores to plot as lines
pns <- clusterMaker(asm_chr, asm_pos, maxGap = 300)
smooth <-  smoother(y = asm_matrix, x = asm_pos, cluster = pns, smoothFunction = locfitByCluster)
asm_matrix_s <- smooth$fitted

# plot the ASM scores across samples

max_y <- max(asm_matrix[!is.na(asm_matrix)])
min_y <- min(asm_matrix[!is.na(asm_matrix)])

i_n <- grep("normal", asm_cols)
i_a <- grep("adenoma", asm_cols)

pdf(""../figures/dame_chr8_ASM_score.pdf", w=10, h=4)

plot(asm_pos, asm_matrix[,1], type='l', col='blue', ylim=c(min_y, max_y), xlab = 'Position on Chromosome 8', ylab = 'ASM Score')

for (i in i_n[-1]) {
  lines(asm_pos, asm_matrix_s[,i], col='blue')
}
for (i in i_a) {
  lines(asm_pos, asm_matrix_s[,i], col='red')
}
for (i in i_n) {
  points(asm_pos, asm_matrix[,i], pch=20, cex=0.5, col='blue')
}
for (i in i_a) {
  points(asm_pos, asm_matrix[,i], pch=20, cex=0.5, col='red')
}

# Shade the DAMR
w <- which(asm_pos >= dame_start & asm_pos <= dame_end)
asm_dame_pos <- asm_pos[w]


polygon(c(asm_dame_pos[1], asm_dame_pos, asm_dame_pos[length(asm_dame_pos)]), 
        c(min_y,rep(max_y, length(asm_dame_pos)),min_y), 
        col=rgb(1, 0, 0,0.05), border=NA)

dev.off()
