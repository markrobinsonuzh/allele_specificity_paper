# allele_specificity_paper
Analysis of allele-specific methylation (and changes in it)

Detailed Methods:

1) Quality Control

Quality control was done on the PE reads. Trimmomatic was subsequently used as follows, with a minimum base quality of 20 for the leading and trailing ends.

$ fastqc sample1_R1.fastq.gz
$ fastqc sample1_R2.fastq.gz

$ java -Xmx2G -cp ~/bin/Trimmomatic-0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 sample1_R1.fastq.gz sample1_R2.fastq.gz sample1_R1_t20l20_paired.fastq.gz sample1_R1_t20l20_single.fastq.gz sample1_R2_t20l20_paired.fastq.gz sample1_R2_t20l20_single.fastq.gz LEADING:20 TRAILING:20 >& sample1_t20l20_trim.out

2) Mapping with Bismark

3) Duplicate Removal

4) Methylation Bias

5) Methtuple


