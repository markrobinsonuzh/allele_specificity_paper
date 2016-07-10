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


pdf(""../figures/DAME_plots.pdf", w=10, h=4)
par(mfrow=c(4,1))

plot_1()
plot_2()



## Functions:

## plot 1: plot ASM scores

plot_1 <- function() {
  
  # load asm_matrix
  load("../data/asm_matrix.rda")
  
  # the rownames of the matrix are called as follows: chr.pos1.pos2
  rows <- rownames(sm)
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
  
  # set x limits for all figures
  x_start <- min(asm_pos)
  x_end <- max(asm_pos)

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
  
  plot(asm_pos, asm_matrix[,1], type='l', col='blue', xlim=c(x_start, x_end), ylim=c(min_y, max_y), xlab = '', ylab = 'ASM Score')
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

}

## plot 2: transformed ASM scores (assumes plot_1 has already run)

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
  
  # Smooth the ASM scores to plot, with lowess by cluster
  pns <- clusterMaker(asm_chr, asm_pos, maxGap = 300)
  smooth <-  smoother(y = asm_t_matrix, x = asm_pos, cluster = pns, smoothFunction = locfitByCluster)
  asm_t_matrix_s <- smooth$fitted
  
  # plot the ASM scores across samples
  max_y <- max(asm_t_matrix[!is.na(asm_t_matrix)])
  min_y <- min(asm_t_matrix[!is.na(asm_t_matrix)])
  
  plot(asm_pos, asm_t_matrix_s[,1], type='l', col='blue', xlim=c(x_start, x_end), ylim=c(max_y,min_y), xlab = '', ylab = 'Transformed ASM Score')
  for (i in i_n[-1]) {
    lines(asm_pos, asm_t_matrix_s[,i], col='blue')
  }
  for (i in i_a) {
    lines(asm_pos, asm_t_matrix_s[,i], col='red')
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
}

## plot 3: plot beta values (t-statistics)





