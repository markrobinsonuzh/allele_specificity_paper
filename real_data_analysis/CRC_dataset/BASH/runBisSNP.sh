#!/bin/bash -x

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/common_all_20170710.vcf.gz

#loop for BisSNP
#.fa file should be sorted and have a .dict file in the same directory

while read filename; do
	name=$(echo "$filename" | egrep -o '[A-Z]{3,4}[1-4]{1}')

	#BisSNP genotyping
	java -Xmx10g -jar BisSNP-1.0.0.jar -R /home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91/GRCh37.91.fa -T BisulfiteGenotyper -I "$filename" -D common_all_20170710.vcf -vfn1 "$name"/"$name".het.snp.raw.vcf -L /home/Shared/steph/130912_HG19_CpGiant_4M_EPI.bed -out_modes EMIT_HET_SNPS_ONLY -nt 20
	 
done < "bam_list.txt"


