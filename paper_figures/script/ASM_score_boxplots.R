# This figures shows the distribution of teh ASM score across the chromosomes in a female and male sample.


# load the list containing the allelicmeth and ASM scores for each sample and combine into one list called 'data'

load("allelicmeth_files_asm_score_data_1.rda")
load("allelicmeth_files_asm_score_data_2.rda")
load("allelicmeth_files_asm_score_data_3.rda")
load("allelicmeth_files_asm_score_data_4.rda")
load("allelicmeth_files_asm_score_data_5.rda")
load("allelicmeth_files_asm_score_data_6.rda")
load("allelicmeth_files_asm_score_data_7.rda")

data <- list()
data <- c(data, data_f_1)
data <- c(data, data_f_2)
data <- c(data, data_f_3)
data <- c(data, data_f_4)
data <- c(data, data_f_5)
data <- c(data, data_f_6)
data <- c(data, data_f_7)

# boxplot ASM distribution

data_plot <- data
for(i in 1:length(data_plot)) {
  s <- data_plot[[i]]
  s$chr <- as.character(s$chr)
  w <- which(s$chr=="chrX" | s$chr=="chrMT" | s$chr=="chrY")
  s[-w,]$chr <- "Autosomal"
  s$chr <- as.factor(s$chr)
  data_plot[[i]] <- s
}

# Female
pdf("./ASM_boxplots_ngs5222.pdf", w=10, h=8)

s <- data_plot[[2]]
boxplot(s$ASM_score ~ s$chr, las=2, col=c(rep("white",2),"orange","white"), outline=F, ylab="ASM Score")

dev.off()

# Male
pdf("./ASM_boxplots_ngs5223.pdf", w=10, h=8)

s <- data_plot[[3]]
boxplot(s$ASM_score ~ s$chr, las=2, col=c(rep("white",2),"orange","white"), outline=F, ylab="ASM Score")

dev.off()


