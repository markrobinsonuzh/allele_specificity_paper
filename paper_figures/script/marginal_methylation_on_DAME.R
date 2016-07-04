# Here we plot the marginal methylation and the smoothed values on a specific DAME on chromosome 8: 141108266-141111068
# Figure 10

# load damr_cov_matrix 
load("methylation_matrix.rda")

library(bumphunter)

damr_chr <- "chr8"
damr_start <- 141108266 
damr_end <- 141111068

# parameters for smoothing
plot_chr <- rep(damr_chr, nrow(damr_cov_matrix))
plot_pos <- as.numeric(limma::strsplit2(rownames(damr_cov_matrix), 'chr8.')[,2])

# sort the positions
o <- order(plot_pos)
plot_pos <- plot_pos[o]
plot_start <- plot_pos[1]
plot_end <- plot_pos[length(plot_pos)]

damr_cov_matrix <- damr_cov_matrix[o,]

plot_pns <- clusterMaker(chr = plot_chr, pos = plot_pos, maxGap = 300)

# use the smoothing function from bumphunter to get smoothed methylation estimates
# this is just to help visualize the marginal methylation per sample

smooth <-  locfitByCluster(damr_cov_matrix, x = plot_pos, cluster = plot_pns)

pdf("./with_bumphunter_package/dame_plots/dame_chr8_marginal_methylation.pdf", w=10, h=8)
plot(plot_pos, smooth$fitted[,1], type='l', col='blue', ylim=c(0,100), lwd=2, 
     main = "Marginal Methylation", xlab = "Position on chr8", ylab = "Percent Methylation")
for (i in c(2,12)) {
  lines(plot_pos, smooth$fitted[,i], col='blue', lwd=2)
}
for (i in c(3,4,5,6,7,8,9,10,11,13)) {
  lines(plot_pos, smooth$fitted[,i], col='red', lwd=2)
}

for (i in c(1,2,12)) {
  points(plot_pos, damr_cov_matrix[,i], pch=20, cex=0.5, col='blue')
}
for (i in c(3,4,5,6,7,8,9,10,11,13)) {
  points(plot_pos, damr_cov_matrix[,i], pch=20, cex=0.5, col='red')
}

# Shade the DAMR
w <- which(plot_pos >= damr_start & plot_pos <= damr_end)
meth_pos <- plot_pos[w]

# df <- as.data.frame(damr_cov_matrix)
# colnames(df) <- limma::strsplit2(colnames(df),"-")[,2]
# df$pos <- plot_pos

# library(reshape2)
# df.m <- melt(df,id.vars='pos', measure.vars=c('normal','Freq.1','Freq.2'))

# p <- ggplot(mry, aes(x=pos, y=, group=rating))
# p + geom_line()


polygon(c(meth_pos[1], meth_pos, meth_pos[length(meth_pos)]), c(0,rep(100, length(meth_pos)),0), 
        col=rgb(1, 0, 0,0.1), border=NA)


dev.off()






