# Pipeline

Raw data can be obtained from ArrayExpress ID: E-MTAB-6949. 

The data was already available form previous analysis, so the downloading of the data is not included here.

Scripts in order of execution:

`BASH/` 
1. `runBismark.sh` align fastq files with bismark. (single job, SLURM script)
2. `runMethtuple.sh` re-sort bam files, and run methtuple. (single job, SLURM script)

`R/`

 3. `fullCancerRun_with2plussamples.R` run DAMEfinder, plot diagnostics.