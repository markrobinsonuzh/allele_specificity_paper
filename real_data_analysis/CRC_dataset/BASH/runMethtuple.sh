#!/bin/sh
#
#SBATCH --job-name=methtuple
#SBATCH -o /home/ubuntu/data/slurm.bismark_%j.out
#SBATCH --ntasks=1

/home/ubuntu/samtools-1.9/samtools sort -n -@ 10 -m 5G -O bam -T $1/_tmp -o $1/$2_qs.bam $1/$2_pe.dedupl_s.bam

methtuple --sc --gzip -m 2 $1/$2_qs.bam 

