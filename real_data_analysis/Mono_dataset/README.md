# Pipeline

Raw data was obtained after the granted permission from the Blueprint Data Access Committee.

Metadata can be obtained from European Genome-phenome Archive (EGA) ID: EGAD00001002523.

Scripts in order of execution:

`BASH/` 
1. `download.sh` to get unaligned crypted bam files.
2. `makeFASTQs.sh` decrypt files and transform to fastq format.
3. `runBismark.sh` align fastq files with bismark. (single job, SLURM script)
4. `mergeBams.sh` merge all bam files from technical replicates. (single job, SLURM script)
5. `runMethtuple.sh` subset the merged bam files by a specified chromosome, and run methtuple. (single job, SLURM script)

`R/`

 6. `fullFemaleMaleRun_methtuple.R` run DAMEfinder, plot diagnostics.