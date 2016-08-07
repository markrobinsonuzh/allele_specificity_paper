# Here we plot the four figures on one of our DAMEs on chromosome 8: 141108266-141111068.
# plot 1: ASM scores
# plot 2: transformed ASM scores
# plot 3: beta(t) plot
# plot 4: marginal methylation
# This is currently figure 10


library(bumphunter)

dame_chr <- "chr8"
dame_start <- 141108266 
dame_end <- 141111068

### load the ASM matrix and Methylation matrix and set the x axis limits for all 4 plots
### The first 3 plots have the same positions but the 4th (methylation plot) has different positions

## ASM Matrix

# load 'asm_matrix' and order it
load("../data/asm_matrix.rda")

# the rownames of the matrix are called as follows: chr.pos1.pos2
rows <- rownames(asm_matrix)
asm_chr <- limma::strsplit2(rows, ".", fixed=T)[,1]
asm_pos1 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,2])
asm_pos2 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,3])

# take the midpoint of pos1 and pos2 as the genomic position for each row
midpt <- floor((asm_pos2 - asm_pos1)/2)
asm_pos <- asm_pos1 + midpt

# make sure the positions are ordered 
o <- order(asm_pos)
asm_matrix <- asm_matrix[o,]
asm_pos <- asm_pos[o]
asm_chr <- asm_chr[o]

## Beta Values

# the cutoff value is obtained using all the loci of the real data set.
# load the beta value vector of the real data set
load("..data/real_beta.rda")

names <- names(beta)
beta_chr <- limma::strsplit2(names, ".", fixed=T)[,1]
beta_pos1 <- as.numeric(limma::strsplit2(names, ".", fixed=T)[,2])
beta_pos2 <- as.numeric(limma::strsplit2(names, ".", fixed=T)[,3])
beta_midpt <- floor((beta_pos2 - beta_pos1)/2)
beta_pos <- beta_pos1 + beta_midpt

# load smoothed beta values
load("../data/real_beta_smoothed.rda")

which_beta <- which(beta_chr==dame_chr & beta_pos >= x_start & beta_pos <= x_end)
beta_pos_p <- beta_pos[which_beta]
beta_chr_p <- beta_chr[which_beta]
beta_p <- beta[which_beta]
smoothed_beta_p <- smoothed_beta[which_beta]

## Methylation Matrix 
# load 'damr_cov_matrix' 
load("../data/methylation_matrix.rda")

meth_chr <- rep(dame_chr, nrow(damr_cov_matrix))
meth_pos <- as.numeric(limma::strsplit2(rownames(damr_cov_matrix), 'chr8.')[,2])

# sort the positions
meth_o <- order(meth_pos)
meeth_pos <- meth_pos[o]
meth_start <- meth_pos[1]
meth_end <- meth_pos[length(meth_pos)]

damr_cov_matrix <- damr_cov_matrix[o,]


## set x limits for all figures
x_start <- min(asm_pos, meth_start)
x_end <- max(asm_pos, meth_end)


### Plot the figures

pdf("../figures/DAME_plots.pdf", w=20, h=16)
par(mfrow=c(4,1))

plot_1()
plot_2()
plot_3()
plot_4()

dev.off()

### Functions:

## plot 1: plots ASM scores

plot_1 <- function() {
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
  
  plot(asm_pos, asm_matrix[,1], type='l', col='blue', xlim=c(x_start, x_end), ylim=c(min_y, max_y), xlab = '', ylab = 'ASM Score', lwd=2, cex.lab=1.3)
  for (i in i_n[-1]) {
    lines(asm_pos, asm_matrix_s[,i], col='blue', lwd = 2)
  }
  for (i in i_a) {
    lines(asm_pos, asm_matrix_s[,i], col='red', lwd = 2)
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
  
  # add legend
  legend(x = x_start, y = 3, c("normal", "adenoma"), 
         lwd = c(2,2), col = c('blue', 'red'), bty = 'n')
  
}

## plot 2: plots transformed ASM scores 
plot_2 <- function() {
  
  # get the ASM_transformed matrix
  modulus_sqrt <- function(values) {
    t_values <- abs(values)
    ret <- sign(values)*sqrt(t_values)
    return(ret)
  }
  
  asm_t_matrix <- matrix(data=NA, nrow=nrow(asm_matrix), ncol=ncol(asm_matrix))
  colnames(asm_t_matrix) <- colnames(asm_matrix)
  rownames(asm_t_matrix) <- rownames(asm_matrix)
  for (i in 1:ncol(asm_t_matrix)) {
    asm_t_matrix[,i] <- modulus_sqrt(asm_matrix[,i])
  }
  
  asm_t_cols <- colnames(asm_t_matrix)
  
  i_n <- grep("normal", asm_t_cols)
  i_a <- grep("adenoma", asm_t_cols)
  
  # Smooth the ASM scores to plot, with lowess by cluster
  pns <- clusterMaker(asm_chr, asm_pos, maxGap = 300)
  smooth <-  smoother(y = asm_t_matrix, x = asm_pos, cluster = pns, smoothFunction = locfitByCluster)
  asm_t_matrix_s <- smooth$fitted
  
  # plot the ASM scores across samples
  max_y <- max(asm_t_matrix[!is.na(asm_t_matrix)])
  min_y <- min(asm_t_matrix[!is.na(asm_t_matrix)])
  
  plot(asm_pos, asm_t_matrix_s[,1], type='l', col='blue', xlim=c(x_start, x_end), ylim=c(min_y,max_y), 
       xlab = '', ylab = 'Transformed ASM Score', lwd=2, cex.lab=1.3)
  for (i in i_n[-1]) {
    lines(asm_pos, asm_t_matrix_s[,i], col='blue', lwd = 2)
  }
  for (i in i_a) {
    lines(asm_pos, asm_t_matrix_s[,i], col='red', lwd = 2)
  }
  for (i in i_n) {
    points(asm_pos, asm_t_matrix[,i], pch=20, cex=0.5, col='blue')
  }
  for (i in i_a) {
    points(asm_pos, asm_t_matrix[,i], pch=20, cex=0.5, col='red')
  }
  
  # Shade the DAMR
  w <- which(asm_pos>=dame_start & asm_pos<=dame_end)
  polygon(c(asm_pos[w][1], asm_pos[w], asm_pos[w][length(asm_pos[w])]), 
          c(min_y,rep(max_y, length(asm_pos[w])),min_y), 
          col=rgb(1, 0, 0,0.05), border=NA)
  # add legend
  legend(x = x_start, y = 1.5, c("normal", "adenoma"), 
         lwd = c(2,2), col = c('blue', 'red'), bty = 'n')
}

## plot 3: plots beta values (t-statistics)
plot_3 <- function() {

  # set cutoff using all beta values
  Q <- 0.9
  cutoff = quantile(abs(smoothed_beta),Q, na.rm=T)
  
  # plot beta and smoothed beta
  y_values <- beta_p[which(!is.na(beta_p))]
  max_y <- max(y_values, cutoff)
  min_y <- min(y_values, -cutoff)
  plot(beta_pos_p, beta_p, xlim = c(x_start, x_end), ylim=c(min_y, max_y), type = 'l', 
       col = 'orange4', xlab = '', ylab = 't-Statistics', lwd = 2, cex.lab=1.3)
  lines(beta_pos_p, smoothed_beta_p, col = 'darkgreen', lwd = 2)
  abline(h = cutoff, col = 'black', lwd = 2)
  abline(h = -cutoff, col = 'black', lwd = 2)
  
  # Shade the DAMR
  which_dame <- which(beta_pos_p>=dame_start & beta_pos_p<=dame_end)
  polygon(c(beta_pos_p[which_dame][1], beta_pos_p[which_dame], beta_pos_p[which_dame][length(beta_pos_p[which_dame])]), 
          c(min_y,rep(max_y, length(beta_pos_p[which_dame])),min_y), 
          col=rgb(1, 0, 0,0.05), border=NA)
  
  # add legend
  legend(x = x_start, y = -50, c("t-Statistic", "Smoothed t-Statistic", "Cutoff"), 
         lwd = c(2,2,2), col = c('orange4', 'darkgreen', 'black'), bty = 'n')
  
}

## plot 4: plots marginal methylation
plot_4 <- function() {
  
  plot_pns <- clusterMaker(chr = meth_chr, pos = meth_pos, maxGap = 300)
  
  # use the smoothing function from bumphunter to get smoothed methylation estimates
  # this is just to help visualize the marginal methylation per sample
  
  smooth <-  smoother(y = damr_cov_matrix, x = meth_pos, cluster = plot_pns, smoothFunction = locfitByCluster)
  
  plot(meth_pos, smooth$fitted[,1], type='l', col='blue', ylim=c(0,100), xlim=c(x_start, x_end), lwd=2, 
       xlab = "Position on chr8", ylab = "Percent Methylation", cex.lab=1.3)
  for (i in c(2,12)) {
    lines(meth_pos, smooth$fitted[,i], col='blue', lwd=2)
  }
  for (i in c(3,4,5,6,7,8,9,10,11,13)) {
    lines(meth_pos, smooth$fitted[,i], col='red', lwd=2)
  }
  
  for (i in c(1,2,12)) {
    points(meth_pos, damr_cov_matrix[,i], pch=20, cex=0.5, col='blue')
  }
  for (i in c(3,4,5,6,7,8,9,10,11,13)) {
    points(meth_pos, damr_cov_matrix[,i], pch=20, cex=0.5, col='red')
  }
  
  # Shade the DAMR
  w <- which(meth_pos >= dame_start & meth_pos <= dame_end)
  meth_pos_p <- meth_pos[w]
  polygon(c(meth_pos_p[1], meth_pos_p, meth_pos_p[length(meth_pos_p)]), c(0,rep(100, length(meth_pos_p)),0), 
          col=rgb(1, 0, 0,0.1), border=NA)
  # add legend
  legend(x = 141106200, y = 40, c("normal", "adenoma"), 
         lwd = c(2,2), col = c('blue', 'red'), bty = 'n')
  
}
