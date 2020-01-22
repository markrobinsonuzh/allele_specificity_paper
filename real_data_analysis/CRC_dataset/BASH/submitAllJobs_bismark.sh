#!/bin/sh 

while read line
do
	files2=$(echo "$line" | sed 's/_R1/_R2/g' | sed 's/val_1/val_2/g')
	#header=$(echo "$line" | cut -d"/" -f6 | sed 's/_R1.fastq.gz//g')
	sample=$(echo "$line" | cut -d"/" -f4 | cut -d"_" -f1-2 ) 
	#header="C001UY"
	sbatch singleJob_bismark.sh "$line" "$files2" "$sample"
done < "meth_bam_files"
