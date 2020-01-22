#!/bin/sh 

while read line 
do
	name=$(echo "$line" | cut -d"/" -f5)
	folder=$(echo "$line" | cut -d"/" -f1-5)
	
	sbatch methtuple_run.sh "$folder" "$name"_single "$line"

done < "ToMethtuple.txt"
