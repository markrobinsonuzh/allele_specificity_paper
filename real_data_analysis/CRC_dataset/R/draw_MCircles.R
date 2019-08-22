## play with MCircle plots
bam_files <- metadata$V3
vcf_files <- metadata$V4 
sample_names <- metadata$V1
reference_file <- "/home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91/GRCh37.91.fa"

# snp2 <- GRanges(9, IRanges(66493966, width = 1))
# snp2 <- GRanges(9, IRanges(66494035, width = 1))
# dame <- GRanges(9,IRanges(66493900,66494060))

snp2 <- GRanges(9, IRanges(99984349, width = 1))
dame <- GRanges(9,IRanges(99984206,99984364))
#dame <- GRanges(9, IRanges(99984177,99984335))

assay(derASM,"snp.table")[subjectHits(findOverlaps(dame, rowRanges(derASM))),]

path <- "/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/"

C2 <- methyl_circle_plot(snp2, vcfFile = gsub("/home/",path, metadata$V4[2]),
                         bamFile = gsub("/home/",path,metadata$V3[2]), 
                         refFile = gsub("/home/",path,reference_file),
                         sampleName = "CRC2",
                         dame = dame,
                         pointSize = 1)

N2 <- methyl_circle_plot(snp2, vcfFile = gsub("/home/",path,metadata$V4[8]),
                         bamFile = gsub("/home/",path,metadata$V3[8]), 
                         refFile = gsub("/home/",path,reference_file),
                         sampleName = "NORM2",
                         dame = dame,
                         pointSize = 1)


C4 <- methyl_circle_plot(snp2, vcfFile = gsub("/home/",path, metadata$V4[4]),
                         bamFile = gsub("/home/",path,metadata$V3[4]), 
                         refFile = gsub("/home/",path,reference_file),
                         sampleName = "CRC4",
                         dame = dame,
                         pointSize = 1)
#doesnt have that SNP ...weird
C4 <- plot.new()

N4 <- methyl_circle_plot(snp2, vcfFile = gsub("/home/",path,metadata$V4[10]),
                         bamFile = gsub("/home/",path,metadata$V3[10]), 
                         refFile = gsub("/home/",path,reference_file),
                         sampleName = "NORM4",
                         dame = dame,
                         pointSize = 1)



C6 <- methyl_circle_plot(snp2, vcfFile = gsub("/home/",path, metadata$V4[6]),
                         bamFile = gsub("/home/",path,metadata$V3[6]), 
                         refFile = gsub("/home/",path,reference_file),
                         sampleName = "CRC6",
                         dame = dame,
                         pointSize = 1)

N6 <- methyl_circle_plot(snp2, vcfFile = gsub("/home/",path,metadata$V4[12]),
                         bamFile = gsub("/home/",path,metadata$V3[12]), 
                         refFile = gsub("/home/",path,reference_file),
                         sampleName = "NORM6",
                         dame = dame,
                         pointSize = 1)


full <- cowplot::plot_grid(C2,C4,C6,N2,N4,N6, ncol=3, nrow = 2, 
                           labels = c("C2", "C4","C6",
                                      "N2","N4","N6"))
  
ggplot2::ggsave("MCircle_plots/MethylcirclesSNP_allCIMPs.png", full, width = 14, height = 12)  
#ggplot2::ggsave("MCircle_plots/MethylcirclesSNP_anotherexample.png", m4, width = 12, height = 10)


### noSNP plots ####
cpgsite <- GRanges(9, IRanges(99984270, width = 1))

m1 <- methyl_circle_plotCpG(cpgsite,
                         bamFile = gsub("/home/",path,metadata$V3[6]), 
                         refFile = gsub("/home/",path,reference_file),
                         dame = dame,
                         pointSize = 1,
                         order = TRUE)

m2 <- methyl_circle_plotCpG(cpgsite,
                            bamFile = gsub("/home/",path,metadata$V3[12]), 
                            refFile = gsub("/home/",path,reference_file),
                            dame = dame,
                            pointSize = 1,
                            order = TRUE)
m4 <- cowplot::plot_grid(m2,m1, ncol=1, nrow = 2, labels = c("NORM","CRC"))
m4
ggplot2::ggsave("MCircle_plots/MethylcirclesCpG_anotherexample.png", m4, width = 12, height = 10)

