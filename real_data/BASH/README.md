## Quality Control

Use the fastqc tool to do quality control and then use the Trimmomatic tool to trim the ends of the reads accordingly.

## Mapping

The BS-seq reads are mapped to the reference genome hg19 using bismark. Duplicate reads are removed with the BLANK command from bismark. Methylationa dn coverage information are also exctracted with the BLANk and BLANK commands.

