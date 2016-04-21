# This figures shows the distribution of teh ASM score across the chromosomes in a female and male sample.


# load the list containing the allelicmeth and ASM scores for each sample and combine into one list called 'data'

load("allelicmeth_files_asm_score_data_1.rda")
load("allelicmeth_files_asm_score_data_2.rda")
load("allelicmeth_files_asm_score_data_3.rda")
load("allelicmeth_files_asm_score_data_4.rda")
load("allelicmeth_files_asm_score_data_5.rda")
load("allelicmeth_files_asm_score_data_6.rda")

data <- list()
data <- c(data, data_f_1)
data <- c(data, data_f_2)
data <- c(data, data_f_3)
data <- c(data, data_f_4)
data <- c(data, data_f_5)
data <- c(data, data_f_6)
