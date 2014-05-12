ngs
===

MokaSubmission.sh
This script supplies the commands such as the number of jobs (samples) and the path names of the output folders to be generated for the MokaPipeline12-05-2014.sh script.


MokaPipeline12-05-2014.sh
This is the latest script for the NGS pipeline. This script is configured to work on the SGE multi-paralell environment and requires commands provided by the MokaSubmission.sh script

Bedfiles folder
Contains:
BedfileInfo.txt relates bedfile to its associated Logfile
BRCA1+2only_LogFile.txt Logfile for how python code was run to generate BRCAdata.bed
BRCAdata.bed bed file for BRCA 1 and 2 panel
BRCA_NGS_CNV_exons_plus_40bpJoCbed.txt BRCA bed file used by Jo Campbell to analyse the BRCA panels
ClinExdata.bed bed file for Clinical exome panel
GSD_all_exons_hg19.bed bed file for the GSD panel
Homo_sapiens.GRCh37.73.dna.primary_assembly.fa Fasta file for Human Index seaquence used for alignment
Homo_sapiens.GRCh37.73.dna.primary_assembly.fa.fai Indexed Fasta file for Human Index seaquence used for alignment
Index.nix ???
PTENdata.bed bed file for PTEN gene panel
PTENonly_LogFile.txt Logfile for how python code was run to generate PTENdata.bed
STK11data.bed bed file for STK11 gene panel
STK11only_LogFile.txt Logfile for how python code was run to generate STK11data.bed
TP53data.bed bed file for TP53 gene panel
TP53only_LogFile.txt Logfile for how python code was run to generate TP53data.bed


