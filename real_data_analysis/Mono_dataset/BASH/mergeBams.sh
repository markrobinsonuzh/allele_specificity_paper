#!/bin/sh
#
#SBATCH --job-name=check_bams
#SBATCH -o /home/ubuntu/data/slurm.bismark_%j.out
#SBATCH --ntasks=5

/home/ubuntu/samtools-1.9/samtools merge -@ 10 $1/$2_merged.bam $3 
/home/ubuntu/samtools-1.9/samtools sort -@ 10 -m 5G -O bam -T _tmp -o $1/$2_merged_s.bam $1/$2_merged.bam
/home/ubuntu/samtools-1.9/samtools index $1/$2_merged_s.bam

