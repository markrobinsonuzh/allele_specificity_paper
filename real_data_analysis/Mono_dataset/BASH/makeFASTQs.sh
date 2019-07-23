#!/bin/bash -x

#Convert unaligned bam files to fastq

while read line; do
	name=$(echo "$line" | cut -f1 -d"." | cut -f3-5 -d"_")
	uncrypt=$(echo "$line" | cut -f1-3 -d".")
	folder=$(echo "$line" | cut -d"/" -f1-2)
	java -jar decryptor.jar EGAdecryptionKey.txt "$line" 
	java -jar /home/Shared_taupo/steph/src/picard/build/libs/picard.jar SamToFastq I="$uncrypt" F="$folder"/"$name"_R1.fastq F2="$folder"/"$name"_R2.fastq 
	gzip "$folder"/"$name"_R1.fastq "$folder"/"$name"_R2.fastq

done < "Ubam.files.txt"
