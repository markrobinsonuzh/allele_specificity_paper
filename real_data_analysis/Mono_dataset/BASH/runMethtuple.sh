#!/bin/sh
#
#SBATCH --job-name=methtuple
#SBATCH -o /home/ubuntu/data/slurm.bismark_%j.out
#SBATCH --ntasks=1

/home/ubuntu/samtools-1.9/samtools view -b $3 3 > $1/$2_3subset.bam
/home/ubuntu/samtools-1.9/samtools sort -n -@ 10 -m 5G -O bam -T $1/_tmp -o $1/$2_3subset_qs.bam $1/$2_3subset.bam

methtuple --sc --gzip -m 2 $1/$2_3subset_qs.bam 

