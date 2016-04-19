## ASM Score and DAMEfinder

We present the Allele-Specific Methylation (ASM) score and the Differential Allele-specific MEthylation finder (DAMEfinder). The ASM score reflects a measure of ASM on a CpG tuple basis using methylation counts of the tuple: MM, MU, UM, and UU, where M stands for methylated and U stands for unmethylated. DAMEfinder uses the ASM score along with the bump hunting approach to detect regions that are Differetnially Allele-specifically MEthylated (DAMEs).

The directories:

* The simulated_data directory contains the simulated BS-seq data that was generated and all the analysis done on this data.
* The real_data directory contains the analysis done on a real adenoma data set.
















change this later

Analysis of allele-specific methylation (and changes in it)

Detailed Methods:

1) Quality Control

Quality control was done on the PE reads. Trimmomatic was subsequently used as follows, with a minimum base quality of 20 for the leading and trailing ends.

$ fastqc sample1_R1.fastq.gz
$ fastqc sample1_R2.fastq.gz

$ java -Xmx2G -cp ~/bin/Trimmomatic-0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 sample1_R1.fastq.gz sample1_R2.fastq.gz sample1_R1_t20l20_paired.fastq.gz sample1_R1_t20l20_single.fastq.gz sample1_R2_t20l20_paired.fastq.gz sample1_R2_t20l20_single.fastq.gz LEADING:20 TRAILING:20 >& sample1_t20l20_trim.out

2) Mapping with Bismark

$ bismark --bowtie2 -p 4 -o $o /home/Shared/data/annotation/_Archive/Human/genome/GRCH37 -1 sample1_R1_t20l20_paired.fastq.gz -2 sample1_R2_t20l20_paired.fastq.gz

3) Duplicate Removal

$  deduplicate_bismark  -p sample1_R1_paired.fastq.gz_bismark_bt2_pe.sam

4) Use methylation extractor to get per CpG coverage and methylation information

$ bismark_methylation_extractor -p --comprehensive sample1_R1_paired.fastq.gz_bismark_bt2_pe.deduplicated.sam

$ bismark2bedGraph --counts -o sample1_R1_paired.fastq.gz_bismark_bt2_pe.deduplicated.bedGraph CpG_context_sample1_R1_paired.fastq.gz_bismark_bt2_pe.deduplicated.txt


4) Methylation Bias

5) Methtuple


