# We use the sm_t matrix which contains the transformed ASM scores and the bump hunting methods to detect DAMEs.

# setwd("/home/Shared_taupo/data/seq/bisulphite_mirco/FASTQ/trimmed_t20l20/bump_hunting/full_bumphunting_real_data/with_bumphunter_package/")

# load sm_t
# load("real_sm_transformed.rda")

library(limma)
library(bumphunter)

############################################################

# The tuples in sm_t are already sorted by the median of the two CpG positions
# The rownames contain the chr, pos1, and pos2 per ASM score in the matrix
rows <- rownames(sm_t)
chr <- limma::strsplit2(rows, ".", fixed=T)[,1]
pos1 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,2])
pos2 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,3])
midpt <- floor((pos2 - pos1)/2)
pos <- pos1 + midpt










