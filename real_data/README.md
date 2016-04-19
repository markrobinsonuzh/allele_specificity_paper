# Real Data Set

The real data set consists of 13 samples of paired end (PE) bisulfite sequencing (BS-seq) reads. 3 of the 13 samples are those of a normal crypt (from 3 separate female individuals). The other 10 are adenoma samples (mixture of male and female individuals). The BS-seq reads were obtained with SureSelect which targets over 350,000 regions on the genome which include CpG rich regions like the CpG islands, shores and shelves, as well as regions known to be differentially methylated in some diseases.


## BASH

This directory contains the commands used to map the BS-seq reads and generate all files that are used to do the analysis in R.

## R

This directory contains the code as implemented in R to calculate the ASM scores and use DAMEfinder to detect the DAMEs.
