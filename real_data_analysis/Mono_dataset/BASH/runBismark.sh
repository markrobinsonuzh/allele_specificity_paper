#!/bin/sh
#
#SBATCH --job-name=bismark_runs
#SBATCH -o /home/ubuntu/data/slurm.bismark_%j.out
#SBATCH --ntasks=16

#Trim
TrimGalore-0.4.5/trim_galore --fastqc -o data/sorjuela/trimmed --paired "$1" "$2" --clip_R1 5 --clip_R2 5 --three_prime_clip_R1 5 --three_prime_clip_R2 5

#Align
Bismark_v0.20.0/bismark --samtools_path /home/ubuntu/samtools-1.9 --path_to_bowtie /home/ubuntu/bowtie2-2.3.4.3 -p 16 -o /home/ubuntu/data/sorjuela/CRC.bismark.bams/ --temp_dir /home/ubuntu/data/sorjuela/CRC.bismark.bams/ --genome /home/ubuntu/data/sorjuela/reference/Bisulfite_Genome.release91 --rg_tag -1 $1 -2 $2 

#Deduplicate
Bismark_v0.20.0/deduplicate_bismark -p --output_dir /home/ubuntu/data/sorjuela/CRC.bismark.bams/ --samtools_path /home/ubuntu/samtools-1.9/ --bam $1

#Sort and index
/home/ubuntu/samtools-1.9/samtools sort -@ 10 -m 5G -O bam -T _tmp -o /home/ubuntu/data/sorjuela/CRC.bismark.bams/$3_pe.dedupl_s.bam $1
/home/ubuntu/samtools-1.9/samtools index /home/ubuntu/data/sorjuela/CRC.bismark.bams/$3_pe.dedupl_s.bam
