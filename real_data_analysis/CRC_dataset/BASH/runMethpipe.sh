#!/bin/bash -x

#amrfinder from the methpipe pipeline

#fix in sherborne
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/
#export LD_LIBRARY_PATH

while read line; do
	name=$(echo "$line" | cut -d"/" -f6 | cut -d"_" -f1)

	/home/sorjuela/methpipe/methpipe-3.4.3/bin/to-mr -o "$name".mr -m bismark -v "$line"

#for some reason the names of the chrom fasta files have to be 1.fa
	/home/sorjuela/methpipe/methpipe-3.4.3/bin/methstates -c /home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91/all_chroms/ -o "$name".epiread -v "$name".mr

#I think theres a bug that reads the mr files and not the epiread files in allelicmeth
	#/home/sorjuela/methpipe/methpipe-3.4.3/bin/allelicmeth -c /home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91/all_chroms/ -o "$name".allelic -v "$name".epiread
	/home/sorjuela/methpipe/methpipe-3.4.3/bin/allelicmeth -v -c /home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91/all_chroms/ -o "$name".allelic_fix "$name".mr
	/home/sorjuela/methpipe/methpipe-3.4.3/bin/amrfinder -o "$name".amr -c /home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91/all_chroms/ "$name".epiread

done < "bam_names"
