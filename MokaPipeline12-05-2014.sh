#!/bin/bash


##Check fastq files are loaded into the PanelNamedata/ProjectNo folder
##Check Mutation reports are loaded into the NextGeneMutationReports folder
##Scripts to integrate Moka pipeline on SGE cluster
##Set the working directory to the current one
#$ -cwd

##Output STDERR and STDOUT to output file place in NGSpipeline_STDOUT_files in Home directory
#$ -j y
#$ -o /home/ryank/NGSpipeline_STDOUT_files
#$ -e /home/ryank/NGSpipeline_STDOUT_files

##Set the shell 
#$ -S /bin/bash

##Set the name of the Job
#$ -N BRCAdataAnalysisbashAllStrandBiasSettingsOn

##Set the Parallel Environment
#$ -pe threaded 8

##Email me when job starts and finishes
#$ -m be
#$ -M kevin.ryan@gsts.com



#########################################################################################
##  READ MAPPING
#########################################################################################
##Generate the variable 'date' with today's date else grep the date from the first run task. This avoids multiple dates being assigned to different tasks if the entire job takes more than one day



##Load modules Novocraft, picardtools, samtools, tabix, vcftools, annovar, bedtools and python

module load novocraft/3.01.02 picard-tools/1.100 samtools/0.1.19 tabix/0.2.6 vcftools/0.1.11 annovar/Oct2013 bedtools/2.17.0 python/2.7_b

#Change directory to where the pre-analysed fastq files are held
cd /home/ryank/$panel

##Assign fastq files to each run task ignoring the name of the program file
sedcommand=$(expr $SGE_TASK_ID + $(expr $SGE_TASK_ID - 1))p
sedcommand1=$(expr $SGE_TASK_ID + $SGE_TASK_ID)p

taskspecificfile=`cat $HOME/${date}run${runstoday}_$project/fastq/filelist.txt | sed -n "$sedcommand"`
taskspecificfile1=`cat $HOME/${date}run${runstoday}_$project/fastq/filelist.txt | sed -n "$sedcommand1"`
taskspecificname=${taskspecificfile1%_L*_R*_*.fastq}
taskspecificname1=${taskspecificname%_*}
#Assign the correct bedfile to the fastq files being analysed
if [[ $taskspecificname1 == *GSD* ]]
then
bedfile='GSD'
echo $bedfile
elif [[ $taskspecificname1 == *BRCA* ]]
then
bedfile='BRCA'
echo $bedfile
elif [[ $taskspecificname1 == *PTEN* ]]
then
bedfile='PTEN'
echo $bedfile
elif [[ $taskspecificname1 == *STK11* ]]
then
bedfile='STK11'
echo $bedfile
elif [[ $taskspecificname1 == *TP53* ]]
then
bedfile='TP53'
echo $bedfile
elif [[ $taskspecificname1 == *ClinEx* ]]
then
bedfile='ClinEx'
echo $bedfile
else
echo "there is no corresponding panel"
fi


echo 'this is the name of the bedfile to be used '${bedfile}data.bed

echo "===================================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "Task index number : $SGE_TASK_ID"
echo "===================================================================="



##Run Moka pipeline


##change directory to working directory for the analysis run
cd $HOME/${date}run${runstoday}_$project


novoalign -c 7 -d $HOME/Databases/Index.nix -f $HOME/${date}run${runstoday}_$project/fastq/$taskspecificfile $HOME/${date}run${runstoday}_$project/fastq/$taskspecificfile1 --Q2Off -F STDFQ -i 200 30 -o SAM -o SoftClip -k -a -g 65 -x 7 2> $HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}_mapping.stats | samtools view -bS - > $HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}.bam



## sort bam file
java -Xmx20g -jar /apps/picardtools/1.100/SortSam.jar INPUT=$HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}.bam OUTPUT=$HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}sorted.bam SORT_ORDER=coordinate TMP_DIR=/scratch/ryank/tmp1 VALIDATION_STRINGENCY=SILENT



## mark duplicates
java -Xmx20g -jar /apps/picardtools/1.100/MarkDuplicates.jar INPUT=$HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}sorted.bam OUTPUT=$HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}dupemarked.bam VALIDATION_STRINGENCY=SILENT METRICS_FILE=$HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}metrics.out TMP_DIR=/scratch/ryank/tmp2 REMOVE_DUPLICATES=TRUE



##Index sort the alignment
samtools index $HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}dupemarked.bam




#########################################################################################
##  VARIANT CALLING
#########################################################################################


##Generate the vcf file for all positions and just the variants
## generates a named pipe called bcfinput 
mkfifo bcfinput$taskspecificname

##Starting with the second command below this note beginning with the "samtools view" command. As STDOUT is generated it is piped to the tee command which simultaneously outputs it to the downstream "bcftools view" command and to the named pipe "bcfinput". This enables STDOUT to be simultaneously processed by the two different bcftools view commands with settings -Nvg and -NG.
bcftools view -Ng bcfinput$taskspecificname - > $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}Total.vcf &

samtools view -bq 20 -F 1796 $HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}dupemarked.bam | samtools mpileup -EDSgu -d 5000 -f $HOME/Databases/Homo_sapiens.GRCh37.73.dna.primary_assembly.fa -L 5000 -F 0.05 -| tee bcfinput$taskspecificname | bcftools view -Nvg - > $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}Variants.vcf



#Generate mpileup file
samtools view -bq 20 -F 1796 $HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}dupemarked.bam | samtools mpileup -EDS -d 5000 -f $HOME/Databases/Homo_sapiens.GRCh37.73.dna.primary_assembly.fa -L 5000 -F 0.05 -> $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}mpileup.txt




##Add sample name to vcf file
sed "s:FORMAT\t-:FORMAT\t${taskspecificname}:g" $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}Variants.vcf > $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}Variants1.vcf


## compress and index
bgzip -f $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}Variants1.vcf
tabix -f $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}Variants1.vcf.gz


## make second vcf just containing coding regions 
vcftools --bed $HOME/Databases/${bedfile}data.bed --gzvcf $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}Variants1.vcf.gz --out $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome --recode --recode-INFO-all


##Filter out variants with a read depth less than 30
vcfutils.pl varFilter -d 30 $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome.recode.vcf > $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome_DP30.vcf


##Filter out variants with a Q value less than 20
vcftools --minQ 20 --vcf $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome_DP30.vcf --recode --recode-INFO-all --out $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome_DP30_minQ20


##compress filtered and unfiltered targetted exome vcf files
bgzip -f $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome.recode.vcf
bgzip -f $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome_DP30.vcf
bgzip -f $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome_DP30_minQ20.recode.vcf


##index filtered and unfiltered targetted exome vcf files
tabix -f $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome.recode.vcf.gz
tabix -f $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome_DP30.vcf.gz
tabix -f $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome_DP30_minQ20.recode.vcf.gz


##decompress gzipped vcf files
bgzip -d $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome_DP30_minQ20.recode.vcf.gz


##############################################################
## Variant reporting
##############################################################


##convert vcf file to annovar format annotate if HOM/HET and PASS/LOW_QUALITY variants
vcfutils.pl varFilter -d 30 $HOME/${date}run${runstoday}_$project/VCF_files/${taskspecificname}targetted_exome_DP30_minQ20.recode.vcf | convert2annovar.pl stdin -format vcf4 -includeinfo | \
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,".",$4,$5,$11,$12,$13,$14,$15}' | \
grep -v "#" | awk 'BEGIN {FS=OFS="\t"} {if (($6 >= 20) && (substr($10,0,3) == "1/1") && (substr($8,1,5) != "INDEL")) print $1,$2,$2,$4,$5,"HOM","PASS",$6,$8,$9,$10 ;
                                else if (($6 < 20) && (substr($10,0,3) == "1/1") && (substr($8,1,5) != "INDEL")) print $1,$2,$2,$4,$5,"HOM","LOW_QUALITY",$6,$8,$9,$10 ;
                                else if (($6 >= 20) && (substr($10,0,3) == "1/1") && (substr($8,1,5) == "INDEL") && ($4 == "-")) print $1,$2,$2,$4,$5,"HOM","PASS",$6,$8,$9,$10 ;
                                else if (($6 < 20) && (substr($10,0,3) == "1/1") && (substr($8,1,5) == "INDEL") && ($4 == "-")) print $1,$2,$2,$4,$5,"HOM","LOW_QUALITY",$6,$8,$9,$10 ;
                                else if (($6 >= 20) && (substr($10,0,3) == "1/1") && (substr($8,1,5) == "INDEL") && ($4 != "-")) print $1,$2,$2+length($4)-1,$4,$5,"HOM","PASS",$6,$8,$9,$10 ;
                                else if (($6 < 20) && (substr($10,0,3) == "1/1") && (substr($8,1,5) == "INDEL") && ($4 != "-")) print $1,$2,$2+length($4)-1,$4,$5,"HOM","LOW_QUALITY",$6,$8,$9,$10 ;
                                else if (($6 >= 20) && (substr($10,0,3) == "0/1") && (substr($8,1,5) != "INDEL")) print $1,$2,$2,$4,$5,"HET","PASS",$6,$8,$9,$10 ;
                                else if (($6 < 20) && (substr($10,0,3) == "0/1") && (substr($8,1,5) != "INDEL")) print $1,$2,$2,$4,$5,"HET","LOW_QUALITY",$6,$8,$9,$10 ;
                                else if (($6 >= 20) && (substr($10,0,3) == "0/1") && (substr($8,1,5) == "INDEL") && ($4 == "-")) print $1,$2,$2,$4,$5,"HET","PASS",$6,$8,$9,$10 ;
                                else if (($6 < 20) && (substr($10,0,3) == "0/1") && (substr($8,1,5) == "INDEL") && ($4 == "-")) print $1,$2,$2,$4,$5,"HET","LOW_QUALITY",$6,$8,$9,$10 ;
                                else if (($6 >= 20) && (substr($10,0,3) == "0/1") && (substr($8,1,5) == "INDEL") && ($4 != "-")) print $1,$2,$2+length($4)-1,$4,$5,"HET","PASS",$6,$8,$9,$10 ;
                                else if (($6 < 20) && (substr($10,0,3) == "0/1") && (substr($8,1,5) == "INDEL") && ($4 != "-")) print $1,$2,$2+length($4)-1,$4,$5,"HET","LOW_QUALITY",$6,$8,$9,$10 ;
        }'  > $HOME/${date}run${runstoday}_$project/Annotation_files/${taskspecificname}.annovar



cd $HOME/${date}run${runstoday}_$project/Annotation_files

## annotate with respect to genes
annotate_variation.pl --buildver hg19 --splicing_threshold 10 -geneanno -exonicsplicing ${taskspecificname}.annovar /apps/annovar/Oct2013/humandb/




## rejig format back to input format
awk 'BEGIN {FS=OFS="\t"} {if ($1 == "exonic;splicing") print $3,$4,$5,$6,$7,$8,substr($1,1,index($1,";")-1),substr($2,1,index($2,";")-1),$9,$10,$11,$12,$13 ;
                     else if (($1 == "splicing" || $1 == "ncRNA_splicing") && index($2,"(") > 0 ) print $3,$4,$5,$6,$7,$8,$1,substr($2,1,index($2,"(")-1),$9,$10,$11,$12,$13 ;
                     else print $3,$4,$5,$6,$7,$8,$1,$2,$9,$10,$11,$12,$13 }' ${taskspecificname}.annovar.variant_function > ${taskspecificname}.annovar.variantsummary




## re-annotate with respect to genes. By rejigging the columns as above you remove duplicated columns
annotate_variation.pl --buildver hg19 --splicing_threshold 10 -geneanno -exonicsplicing ${taskspecificname}.annovar.variantsummary /apps/annovar/Oct2013/humandb/



##Rejig format back to input
awk 'BEGIN {FS=OFS="\t"} {if ((substr($1,0,6) != "exonic") && ($1 != "splicing")) print $3,$4,$5,$6,$7,$8,$9,$10,$9,".",$11,$12,$13,$14,$15;      
                     else if ($1 == "splicing") print $3,$4,$5,$6,$7,$8,$9,$10,$9,$10":"substr($2,index($2,"(")+1,index($2,")")-index($2,"(")-1),$11,$12,$13,$14,$15}' ${taskspecificname}.annovar.variantsummary.variant_function > ${taskspecificname}.annovar.variantsummary1
awk 'BEGIN {FS=OFS="\t"} {print $4,$5,$6,$7,$8,$9,$10,$11,$2,$3,$12,$13,$14,$15,$16}' ${taskspecificname}.annovar.variantsummary.exonic_variant_function >>  ${taskspecificname}.annovar.variantsummary1



## cross reference with dbSNP137a
annotate_variation.pl --buildver hg19 --filter -dbtype snp137 ${taskspecificname}.annovar.variantsummary1 /apps/annovar/Oct2013/humandb/



## rejig format back to input format
awk 'BEGIN {FS=OFS="\t"} {print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$2,$13,$14,$15,$16,$17}' ${taskspecificname}.annovar.variantsummary1.hg19_snp137_dropped > ${taskspecificname}.annovar.variantsummary2
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,".",$11,$12,$13,$14,$15}' ${taskspecificname}.annovar.variantsummary1.hg19_snp137_filtered >> ${taskspecificname}.annovar.variantsummary2




## cross reference with EVS
annotate_variation.pl --buildver hg19 --filter -dbtype esp6500si_all ${taskspecificname}.annovar.variantsummary2 /apps/annovar/Oct2013/humandb/ESP




## rejig format back to input format
awk 'BEGIN {FS=OFS="\t"} {print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$2,$14,$15,$16,$17,$18}' ${taskspecificname}.annovar.variantsummary2.hg19_esp6500si_all_dropped > ${taskspecificname}.annovar.variantsummary3 
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,".",$12,$13,$14,$15,$16}' ${taskspecificname}.annovar.variantsummary2.hg19_esp6500si_all_filtered >> ${taskspecificname}.annovar.variantsummary3 



## cross reference with 1000g
annotate_variation.pl --buildver hg19 --filter -dbtype 1000g2012apr_all ${taskspecificname}.annovar.variantsummary3 /apps/annovar/Oct2013/humandb/1000g



## rejig format back to input format
awk 'BEGIN {FS=OFS="\t"} {print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$2,$15,$16,$17,$18,$19}' ${taskspecificname}.annovar.variantsummary3.hg19_ALL.sites.2012_04_dropped > ${taskspecificname}.annovar.variantsummary4
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,".",$13,$14,$15,$16,$17}' ${taskspecificname}.annovar.variantsummary3.hg19_ALL.sites.2012_04_filtered >> ${taskspecificname}.annovar.variantsummary4

##Cross reference with SIFT
annotate_variation.pl --buildver hg19 --filter -dbtype avsift ${taskspecificname}.annovar.variantsummary4 /apps/annovar/Oct2013/humandb/SIFT


## rejig format back to input format
awk 'BEGIN {FS=OFS="\t"} {print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$2,$16,$17,$18,$19,$20}' ${taskspecificname}.annovar.variantsummary4.hg19_avsift_dropped > ${taskspecificname}.annovar.variantsummary5
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,".",$14,$15,$16,$17,$18}' ${taskspecificname}.annovar.variantsummary4.hg19_avsift_filtered >> ${taskspecificname}.annovar.variantsummary5



##assess and report novelty
awk 'BEGIN {FS=OFS="\t"} { if (($10 == ".") && ($11 == ".") && ($12 == ".") && ($13 == ".") && ($14 == ".")) print $1,$2,$3,$4,$5,$6,$7,$8,$9,"NOVEL",$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;
                else if (($10 != ".") || ($11 != ".") || ($12 != ".") || ($13 != ".") || ($14 != ".")) print $1,$2,$3,$4,$5,$6,$7,$8,$9,".",$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; 
        }' ${taskspecificname}.annovar.variantsummary5 > ${taskspecificname}.annovar.variantsummary6



##grep all variants assigned a PASS that correspond to exonic|splicing regions and write to ${taskspecificname}_total_uniq.var 
grep -v LOW_QUALITY $HOME/${date}run${runstoday}_$project/Annotation_files/${taskspecificname}.annovar.variantsummary6 | grep -v ncRNA | egrep 'exonic|splicing' > ${taskspecificname}_total_uniq.var 
awk 'BEGIN {FS=OFS="\t"}{if ($6 == "HET") print $0}' ${taskspecificname}_total_uniq.var > ${taskspecificname}_het_uniq.var 
awk 'BEGIN {FS=OFS="\t"}{if ($6 == "HOM") print $0}' ${taskspecificname}_total_uniq.var > ${taskspecificname}_hom_uniq.var 



## concatanate the file and pipe into command wc -l which counts the number of corresponding lines
total_var_all=`cat ${taskspecificname}_total_uniq.var | wc -l`
total_var_dbSNP=`awk 'BEGIN {FS="\t"} {if ($12 != ".") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_var_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ($12 == ".") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_var_known=`awk 'BEGIN {FS="\t"} {if ($10 != "NOVEL") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_var_novel=`awk 'BEGIN {FS="\t"} {if ($10 == "NOVEL") print $0}' ${taskspecificname}_total_uniq.var | wc -l`

total_coding_var_all=`awk 'BEGIN {FS="\t"} {if ($7 == "exonic") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_var_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_var_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($12 == ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_var_known=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_var_novel=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`


total_coding_SNV_all=`awk 'BEGIN {FS="\t"} {if (substr($9,length($9)-2) == "SNV") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_SNV_known=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_SNV_novel=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`


total_coding_synonymous_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "synonymous SNV") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_synonymous_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_synonymous_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_synonymous_SNV_known=`awk 'BEGIN {FS="\t"} {if (($10 == "synonymous SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_synonymous_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($10 == "synonymous SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`

total_coding_nonsynonymous_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "nonsynonymous SNV") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_nonsynonymous_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_nonsynonymous_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_nonsynonymous_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_nonsynonymous_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`


total_coding_stopgain_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "stopgain SNV") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_stopgain_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_stopgain_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_stopgain_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_stopgain_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`

total_coding_stoploss_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "stoploss SNV") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_stoploss_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_stoploss_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_stoploss_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_stoploss_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`

total_coding_DELETION_all=`awk 'BEGIN {FS="\t"} {if (substr($9,length($9)-7) == "deletion") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_DELETION_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_DELETION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_DELETION_known=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_DELETION_novel=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`


total_coding_INSERTION_all=`awk 'BEGIN {FS="\t"} {if (substr($9,length($9)-8) == "insertion") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_INSERTION_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_INSERTION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_INSERTION_known=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_INSERTION_novel=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`

total_coding_FRAMESHIFT_DELETION_all=`awk 'BEGIN {FS="\t"} {if ($9 == "frameshift deletion") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_FRAMESHIFT_DELETION_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_FRAMESHIFT_DELETION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_FRAMESHIFT_DELETION_known=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_FRAMESHIFT_DELETION_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`

total_coding_FRAMESHIFT_INSERTION_all=`awk 'BEGIN {FS="\t"} {if ($9 == "frameshift insertion") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_FRAMESHIFT_INSERTION_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_FRAMESHIFT_INSERTION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_FRAMESHIFT_INSERTION_known=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_coding_FRAMESHIFT_INSERTION_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`

total_splice_all=`awk 'BEGIN {FS="\t"} {if ($7 == "splicing") print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_splice_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_splice_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_splice_known=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_splice_novel=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`

total_transitions_all=`awk 'BEGIN {FS="\t"} {if ((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_transitions_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_transitions_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($12 == "."))print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_transitions_known=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_transitions_novel=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`


total_transversions_all=`awk 'BEGIN {FS="\t"} {if ((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_transversions_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($12 != ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_transversions_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($12 == ".")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_transversions_known=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($10 != "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`
total_transversions_novel=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($10 == "NOVEL")) print $0}' ${taskspecificname}_total_uniq.var | wc -l`



## scale=2 is the number of digits after the decimal point 
## bc is a language that supports arbitrary precision numbers. A standard maths library is available by command line option

total_ts_tv_ratio_all=$(echo "scale=2; $total_transitions_all/$total_transversions_all" | bc)
total_ts_tv_ratio_dbSNP=$(echo "scale=2; $total_transitions_dbSNP/$total_transversions_dbSNP" | bc)
total_ts_tv_ratio_not_dbSNP=$(echo "scale=2; $total_transitions_not_dbSNP/$total_transversions_not_dbSNP" | bc)
total_ts_tv_ratio_known=$(echo "scale=2; $total_transitions_known/$total_transversions_known" | bc)
total_ts_tv_ratio_novel=$(echo "scale=2; $total_transitions_novel/$total_transversions_novel" | bc)

het_var_all=`cat ${taskspecificname}_het_uniq.var | wc -l`
het_var_dbSNP=`awk 'BEGIN {FS="\t"} {if ($12 != ".") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_var_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ($12 == ".") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_var_known=`awk 'BEGIN {FS="\t"} {if ($10 != "NOVEL") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_var_novel=`awk 'BEGIN {FS="\t"} {if ($10 == "NOVEL") print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_var_all=`awk 'BEGIN {FS="\t"} {if ($7 == "exonic") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_var_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_var_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($12 == ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_var_known=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_var_novel=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_SNV_all=`awk 'BEGIN {FS="\t"} {if (substr($9,length($9)-2) == "SNV") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_SNV_known=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_SNV_novel=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_synonymous_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "synonymous SNV") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_synonymous_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_synonymous_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_synonymous_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_synonymous_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_nonsynonymous_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "nonsynonymous SNV") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_nonsynonymous_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_nonsynonymous_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_nonsynonymous_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_nonsynonymous_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_stopgain_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "stopgain SNV") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_stopgain_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_stopgain_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_stopgain_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_stopgain_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_stoploss_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "stoploss SNV") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_stoploss_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_stoploss_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_stoploss_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_stoploss_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_DELETION_all=`awk 'BEGIN {FS="\t"} {if (substr($9,length($9)-7) == "deletion") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_DELETION_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_DELETION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_DELETION_known=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_DELETION_novel=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_INSERTION_all=`awk 'BEGIN {FS="\t"} {if (substr($9,length($9)-8) == "insertion") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_INSERTION_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_INSERTION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_INSERTION_known=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_INSERTION_novel=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_FRAMESHIFT_DELETION_all=`awk 'BEGIN {FS="\t"} {if ($9 == "frameshift deletion") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_FRAMESHIFT_DELETION_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_FRAMESHIFT_DELETION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_FRAMESHIFT_DELETION_known=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_FRAMESHIFT_DELETION_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_coding_FRAMESHIFT_INSERTION_all=`awk 'BEGIN {FS="\t"} {if ($9 == "frameshift insertion") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_FRAMESHIFT_INSERTION_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_FRAMESHIFT_INSERTION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_FRAMESHIFT_INSERTION_known=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_coding_FRAMESHIFT_INSERTION_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_splice_all=`awk 'BEGIN {FS="\t"} {if ($7 == "splicing") print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_splice_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_splice_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_splice_known=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_splice_novel=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_transitions_all=`awk 'BEGIN {FS="\t"} {if ((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_transitions_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_transitions_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($12 == "."))print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_transitions_known=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_transitions_novel=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_transversions_all=`awk 'BEGIN {FS="\t"} {if ((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_transversions_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($12 != ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_transversions_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($12 == ".")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_transversions_known=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($10 != "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`
het_transversions_novel=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($10 == "NOVEL")) print $0}' ${taskspecificname}_het_uniq.var | wc -l`

het_ts_tv_ratio_all=$(echo "scale=2; $het_transitions_all/$het_transversions_all" | bc)
het_ts_tv_ratio_dbSNP=$(echo "scale=2; $het_transitions_dbSNP/$het_transversions_dbSNP" | bc)
het_ts_tv_ratio_not_dbSNP=$(echo "scale=2; $het_transitions_not_dbSNP/$het_transversions_not_dbSNP" | bc)
het_ts_tv_ratio_known=$(echo "scale=2; $het_transitions_known/$het_transversions_known" | bc)
het_ts_tv_ratio_novel=$(echo "scale=2; $het_transitions_novel/$het_transversions_novel" | bc)


hom_var_all=`cat ${taskspecificname}_hom_uniq.var | wc -l`
hom_var_dbSNP=`awk 'BEGIN {FS="\t"} {if ($12 != ".") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_var_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ($12 == ".") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_var_known=`awk 'BEGIN {FS="\t"} {if ($10 != "NOVEL") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_var_novel=`awk 'BEGIN {FS="\t"} {if ($10 == "NOVEL") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_var_all=`awk 'BEGIN {FS="\t"} {if ($7 == "exonic") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_var_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_var_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($12 == ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_var_known=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_var_novel=`awk 'BEGIN {FS="\t"} {if (($7 == "exonic") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_SNV_all=`awk 'BEGIN {FS="\t"} {if (substr($9,length($9)-2) == "SNV") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_SNV_known=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_SNV_novel=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-2) == "SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_synonymous_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "synonymous SNV") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_synonymous_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_synonymous_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_synonymous_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_synonymous_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "synonymous SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_nonsynonymous_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "nonsynonymous SNV") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_nonsynonymous_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_nonsynonymous_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_nonsynonymous_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_nonsynonymous_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "nonsynonymous SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_stopgain_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "stopgain SNV") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_stopgain_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_stopgain_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_stopgain_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_stopgain_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "stopgain SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_stoploss_SNV_all=`awk 'BEGIN {FS="\t"} {if ($9 == "stoploss SNV") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_stoploss_SNV_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_stoploss_SNV_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_stoploss_SNV_known=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_stoploss_SNV_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "stoploss SNV") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_DELETION_all=`awk 'BEGIN {FS="\t"} {if (substr($9,length($9)-7) == "deletion") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_DELETION_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_DELETION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_DELETION_known=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_DELETION_novel=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-7) == "deletion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_INSERTION_all=`awk 'BEGIN {FS="\t"} {if (substr($9,length($9)-8) == "insertion") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_INSERTION_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_INSERTION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_INSERTION_known=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_INSERTION_novel=`awk 'BEGIN {FS="\t"} {if ((substr($9,length($9)-8) == "insertion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_FRAMESHIFT_DELETION_all=`awk 'BEGIN {FS="\t"} {if ($9 == "frameshift deletion") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_FRAMESHIFT_DELETION_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_FRAMESHIFT_DELETION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_FRAMESHIFT_DELETION_known=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_FRAMESHIFT_DELETION_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift deletion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_coding_FRAMESHIFT_INSERTION_all=`awk 'BEGIN {FS="\t"} {if ($9 == "frameshift insertion") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_FRAMESHIFT_INSERTION_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_FRAMESHIFT_INSERTION_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_FRAMESHIFT_INSERTION_known=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_coding_FRAMESHIFT_INSERTION_novel=`awk 'BEGIN {FS="\t"} {if (($9 == "frameshift insertion") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_splice_all=`awk 'BEGIN {FS="\t"} {if ($7 == "splicing") print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_splice_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_splice_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_splice_known=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_splice_novel=`awk 'BEGIN {FS="\t"} {if (($7 == "splicing") && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_transitions_all=`awk 'BEGIN {FS="\t"} {if ((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_transitions_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_transitions_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($12 == "."))print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_transitions_known=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_transitions_novel=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "G")) || (($4 == "G") && ($5 == "A")) || (($4 == "C") && ($5 == "T")) || (($4 == "T") && ($5 == "C"))) && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_transversions_all=`awk 'BEGIN {FS="\t"} {if ((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_transversions_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($12 != ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_transversions_not_dbSNP=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($12 == ".")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_transversions_known=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($10 != "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`
hom_transversions_novel=`awk 'BEGIN {FS="\t"} {if (((($4 == "A") && ($5 == "C")) || (($4 == "C") && ($5 == "A")) || (($4 == "A") && ($5 == "T")) || (($4 == "T") && ($5 == "A")) || (($4 == "C") && ($5 == "G")) || (($4 == "G") && ($5 == "C")) || (($4 == "G") && ($5 == "T")) || (($4 == "T") && ($5 == "G"))) && ($10 == "NOVEL")) print $0}' ${taskspecificname}_hom_uniq.var | wc -l`

hom_ts_tv_ratio_all=$(echo "scale=2; $hom_transitions_all/$hom_transversions_all" | bc)
hom_ts_tv_ratio_dbSNP=$(echo "scale=2; $hom_transitions_dbSNP/$hom_transversions_dbSNP" | bc)
hom_ts_tv_ratio_not_dbSNP=$(echo "scale=2; $hom_transitions_not_dbSNP/$hom_transversions_not_dbSNP" | bc)
hom_ts_tv_ratio_known=$(echo "scale=2; $hom_transitions_known/$hom_transversions_known" | bc)
hom_ts_tv_ratio_novel=$(echo "scale=2; $hom_transitions_novel/$hom_transversions_novel" | bc)

##print statements for each $variable defined above wrtten to file
printf "$1\t$1\t$1\t$1\t$1\t$1\n" > ${taskspecificname}_annovar_variant.stats
printf "variant_type\tall\tin_dbSNP135\tnot_in_dbSNP135\tknown\tnovel\n" >> ${taskspecificname}_annovar_variant.stats

printf "variants\t$total_var_all\t$total_var_dbSNP\t$total_var_not_dbSNP\t$total_var_known\t$total_var_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_variants\t$het_var_all\t$het_var_dbSNP\t$het_var_not_dbSNP\t$het_var_known\t$het_var_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_variants\t$hom_var_all\t$hom_var_dbSNP\t$hom_var_not_dbSNP\t$hom_var_known\t$hom_var_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "coding_variants\t$total_coding_var_all\t$total_coding_var_dbSNP\t$total_coding_var_not_dbSNP\t$total_coding_var_known\t$total_coding_var_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_coding_variants\t$het_coding_var_all\t$het_coding_var_dbSNP\t$het_coding_var_not_dbSNP\t$het_coding_var_known\t$het_coding_var_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_coding_variants\t$hom_coding_var_all\t$hom_coding_var_dbSNP\t$hom_coding_var_not_dbSNP\t$hom_coding_var_known\t$hom_coding_var_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "splice_variants\t$total_splice_all\t$total_splice_dbSNP\t$total_splice_not_dbSNP\t$total_splice_known\t$total_splice_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_splice_variants\t$het_splice_all\t$het_splice_dbSNP\t$het_splice_not_dbSNP\t$het_splice_known\t$het_splice_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_splice_variants\t$hom_splice_all\t$hom_splice_dbSNP\t$hom_splice_not_dbSNP\t$hom_splice_known\t$hom_splice_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "nonsynonymous_SNVs\t$total_coding_nonsynonymous_SNV_all\t$total_coding_nonsynonymous_SNV_dbSNP\t$total_coding_nonsynonymous_SNV_not_dbSNP\t$total_coding_nonsynonymous_SNV_known\t$total_coding_nonsynonymous_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_nonsynonymous_SNVs\t$het_coding_nonsynonymous_SNV_all\t$het_coding_nonsynonymous_SNV_dbSNP\t$het_coding_nonsynonymous_SNV_not_dbSNP\t$het_coding_nonsynonymous_SNV_known\t$het_coding_nonsynonymous_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_nonsynonymous_SNVs\t$hom_coding_nonsynonymous_SNV_all\t$hom_coding_nonsynonymous_SNV_dbSNP\t$hom_coding_nonsynonymous_SNV_not_dbSNP\t$hom_coding_nonsynonymous_SNV_known\t$hom_coding_nonsynonymous_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "synonymous_SNVs\t$total_coding_synonymous_SNV_all\t$total_coding_synonymous_SNV_dbSNP\t$total_coding_synonymous_SNV_not_dbSNP\t$total_coding_synonymous_SNV_known\t$total_coding_synonymous_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_synonymous_SNVs\t$het_coding_synonymous_SNV_all\t$het_coding_synonymous_SNV_dbSNP\t$het_coding_synonymous_SNV_not_dbSNP\t$het_coding_synonymous_SNV_known\t$het_coding_synonymous_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_synonymous_SNVs\t$hom_coding_synonymous_SNV_all\t$hom_coding_synonymous_SNV_dbSNP\t$hom_coding_synonymous_SNV_not_dbSNP\t$hom_coding_synonymous_SNV_known\t$hom_coding_synonymous_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "stoploss_SNVs\t$total_coding_stoploss_SNV_all\t$total_coding_stoploss_SNV_dbSNP\t$total_coding_stoploss_SNV_not_dbSNP\t$total_coding_stoploss_SNV_known\t$total_coding_stoploss_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_stoploss_SNVs\t$het_coding_stoploss_SNV_all\t$het_coding_stoploss_SNV_dbSNP\t$het_coding_stoploss_SNV_not_dbSNP\t$het_coding_stoploss_SNV_known\t$het_coding_stoploss_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_stoploss_SNVs\t$hom_coding_stoploss_SNV_all\t$hom_coding_stoploss_SNV_dbSNP\t$hom_coding_stoploss_SNV_not_dbSNP\t$hom_coding_stoploss_SNV_known\t$hom_coding_stoploss_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "stopgain_SNVs\t$total_coding_stopgain_SNV_all\t$total_coding_stopgain_SNV_dbSNP\t$total_coding_stopgain_SNV_not_dbSNP\t$total_coding_stopgain_SNV_known\t$total_coding_stopgain_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_stopgain_SNVs\t$het_coding_stopgain_SNV_all\t$het_coding_stopgain_SNV_dbSNP\t$het_coding_stopgain_SNV_not_dbSNP\t$het_coding_stopgain_SNV_known\t$het_coding_stopgain_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_stopgain_SNVs\t$hom_coding_stopgain_SNV_all\t$hom_coding_stopgain_SNV_dbSNP\t$hom_coding_stopgain_SNV_not_dbSNP\t$hom_coding_stopgain_SNV_known\t$hom_coding_stopgain_SNV_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "deletions\t$total_coding_DELETION_all\t$total_coding_DELETION_dbSNP\t$total_coding_DELETION_not_dbSNP\t$total_coding_DELETION_known\t$total_coding_DELETION_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_deletions\t$het_coding_DELETION_all\t$het_coding_DELETION_dbSNP\t$het_coding_DELETION_not_dbSNP\t$het_coding_DELETION_known\t$het_coding_DELETION_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_deletions\t$hom_coding_DELETION_all\t$hom_coding_DELETION_dbSNP\t$hom_coding_DELETION_not_dbSNP\t$hom_coding_DELETION_known\t$hom_coding_DELETION_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "insertions\t$total_coding_INSERTION_all\t$total_coding_INSERTION_dbSNP\t$total_coding_INSERTION_not_dbSNP\t$total_coding_INSERTION_known\t$total_coding_INSERTION_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_insertions\t$het_coding_INSERTION_all\t$het_coding_INSERTION_dbSNP\t$het_coding_INSERTION_not_dbSNP\t$het_coding_INSERTION_known\t$het_coding_INSERTION_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_insertions\t$hom_coding_INSERTION_all\t$hom_coding_INSERTION_dbSNP\t$hom_coding_INSERTION_not_dbSNP\t$hom_coding_INSERTION_known\t$hom_coding_INSERTION_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "frameshift_deletions\t$total_coding_FRAMESHIFT_DELETION_all\t$total_coding_FRAMESHIFT_DELETION_dbSNP\t$total_coding_FRAMESHIFT_DELETION_not_dbSNP\t$total_coding_FRAMESHIFT_DELETION_known\t$total_coding_FRAMESHIFT_DELETION_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_frameshift_deletions\t$het_coding_FRAMESHIFT_DELETION_all\t$het_coding_FRAMESHIFT_DELETION_dbSNP\t$het_coding_FRAMESHIFT_DELETION_not_dbSNP\t$het_coding_FRAMESHIFT_DELETION_known\t$het_coding_FRAMESHIFT_DELETION_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_frameshift_deletions\t$hom_coding_FRAMESHIFT_DELETION_all\t$hom_coding_FRAMESHIFT_DELETION_dbSNP\t$hom_coding_FRAMESHIFT_DELETION_not_dbSNP\t$hom_coding_FRAMESHIFT_DELETION_known\t$hom_coding_FRAMESHIFT_DELETION_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "frameshift_insertions\t$total_coding_FRAMESHIFT_INSERTION_all\t$total_coding_FRAMESHIFT_INSERTION_dbSNP\t$total_coding_FRAMESHIFT_INSERTION_not_dbSNP\t$total_coding_FRAMESHIFT_INSERTION_known\t$total_coding_FRAMESHIFT_INSERTION_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_frameshift_insertions\t$het_coding_FRAMESHIFT_INSERTION_all\t$het_coding_FRAMESHIFT_INSERTION_dbSNP\t$het_coding_FRAMESHIFT_INSERTION_not_dbSNP\t$het_coding_FRAMESHIFT_INSERTION_known\t$het_coding_FRAMESHIFT_INSERTION_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_frameshift_insertions\t$hom_coding_FRAMESHIFT_INSERTION_all\t$hom_coding_FRAMESHIFT_INSERTION_dbSNP\t$hom_coding_FRAMESHIFT_INSERTION_not_dbSNP\t$hom_coding_FRAMESHIFT_INSERTION_known\t$hom_coding_FRAMESHIFT_INSERTION_novel\n" >> ${taskspecificname}_annovar_variant.stats

printf "ts_tv_ratio\t$total_ts_tv_ratio_all\t$total_ts_tv_ratio_dbSNP\t$total_ts_tv_ratio_not_dbSNP\t$total_ts_tv_ratio_known\t$total_ts_tv_ratio_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "het_ts_tv_ratio\t$het_ts_tv_ratio_all\t$het_ts_tv_ratio_dbSNP\t$het_ts_tv_ratio_not_dbSNP\t$het_ts_tv_ratio_known\t$het_ts_tv_ratio_novel\n" >> ${taskspecificname}_annovar_variant.stats
printf "hom_ts_tv_ratio\t$hom_ts_tv_ratio_all\t$hom_ts_tv_ratio_dbSNP\t$hom_ts_tv_ratio_not_dbSNP\t$hom_ts_tv_ratio_known\t$hom_ts_tv_ratio_novel\n" >> ${taskspecificname}_annovar_variant.stats

##########Validation step

# coverage calculations

samtools view -bq 20 -F 1796 $HOME/${date}run${runstoday}_$project/Alignment_Output/${taskspecificname}dupemarked.bam | bamToBed -i stdin > ${taskspecificname}final.bed 

#Shows a break down of the read depth of each specific region in the bed file provided by -b option  

coverageBed -hist -a ${taskspecificname}final.bed -b $HOME/Databases/${bedfile}data.bed > ${taskspecificname}coverage_gencode_interval.bed

#Coverage per exon ie the number of reads per exon
coverageBed -a ${taskspecificname}final.bed -b $HOME/Databases/${bedfile}data.bed | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4,$5}' > ${taskspecificname}coverage_targets.bed


#Calculates read depth per base

coverageBed -d -a ${taskspecificname}final.bed -b $HOME/Databases/${bedfile}data.bed > ${taskspecificname}coverage_targets1.bed 


#efficiency of capture
#Calculate the total number of reads from bed file generated from BAM file
# calculate the total number of tab-separated rows (ie the number of records) for final.bed file and assign it to the variable 'reads'
reads=`awk 'END {OFS = "\t"; print NR}' ${taskspecificname}final.bed`

# calculate the sum of all tab-separated elements in column 4 from coverage_targets.bed file which calculates the total number of reads mapping to the Region of interest
mapped_to_target_reads=`awk '{SUM += $5} END {OFS = "\t";print SUM}' ${taskspecificname}coverage_targets.bed`

#Percentage of reads mapping to the region of interest
# calculate the percentage of mapped reads in coverage_targets.bed file to two decimal places assign the percentage value to the variable percent1

percent1=`awk 'BEGIN{printf("%0.2f", ('$mapped_to_target_reads' / '$reads') * 100)}'`

#coverage gives the range of read depths across the capture region
grep all ${taskspecificname}coverage_gencode_interval.bed > ${taskspecificname}coverage.hist
#Average number of reads overall
meancov=`awk '{if ($2>=1) (SUM += $2*$5)} END {printf ("%0.2f", SUM)}' ${taskspecificname}coverage.hist`

#completeness of coverage (ie how much of target region has been covered to a particular depth)
cov1xpc=`awk '{if ($2>=1) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' ${taskspecificname}coverage.hist`
cov10xpc=`awk '{if ($2>=10) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' ${taskspecificname}coverage.hist`
cov30xpc=`awk '{if ($2>=30) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' ${taskspecificname}coverage.hist`
cov100xpc=`awk '{if ($2>=100) (SUM += $5)} END {printf ("%0.2f", SUM*100)}' ${taskspecificname}coverage.hist`

#Total number of bases covered to a particular level
cov1x=`awk '{if ($2>=1) (SUM += $3)} END {print SUM}' ${taskspecificname}coverage.hist`
cov10x=`awk '{if ($2>=10) (SUM += $3)} END {print SUM}' ${taskspecificname}coverage.hist`
cov30x=`awk '{if ($2>=30) (SUM += $3)} END {print SUM}' ${taskspecificname}coverage.hist`
cov100x=`awk '{if ($2>=100) (SUM += $3)} END {print SUM}' ${taskspecificname}coverage.hist`

#report generation
printf "total_reads\t"$reads"\n" > ${taskspecificname}coverage.stats
printf "mapped_to_target_reads\t"$mapped_to_target_reads"\n" >> ${taskspecificname}coverage.stats
printf "percentage\t"$percent1"\n" >> ${taskspecificname}coverage.stats
printf "mean_coverage\t"$meancov"\n" >> ${taskspecificname}coverage.stats
printf "No of accessible_target_bases at 1x coverage\t"$cov1x"\n" >> ${taskspecificname}coverage.stats
printf "No of accessible_target_bases at 10x coverage\t"$cov10x"\n" >> ${taskspecificname}coverage.stats
printf "Percentage equivalent\t"$cov10xpc"\n" >> ${taskspecificname}coverage.stats
printf "No of accessible_target_bases at 30x coverage\t"$cov30x"\n" >> ${taskspecificname}coverage.stats
printf "Percentage equivalent\t"$cov30xpc"\n" >> ${taskspecificname}coverage.stats
printf "No of accessible_target_bases at 100x coverage\t"$cov100x"\n" >> ${taskspecificname}coverage.stats
printf "Percentage equivalent\t"$cov100xpc"\n" >> ${taskspecificname}coverage.stats

#Don't need the python script for BRCA cross-referencing
#Execute python script to cross reference data generated by Moka with variants generated by NextGene software which went on to be Sanger sequenced
#python $HOME/Sangercheckv1 -m $HOME/${date}run${runstoday}_$project/Annotation_files/${taskspecificname}.annovar.variantsummary -s /home/ryank/NextGeneMutationReports/${taskspecificname1}_MutationReport.csv -o ${taskspecificname}SangerCrossRef.csv

#cd $HOME/NextGeneMutationReports/
#Move the analysed Mutation report to the Annotation folder of the run directory
#mv ${taskspecificname1}_MutationReport.csv $HOME/${date}run${runstoday}_$project/Annotation_files/

#Copy STDOUT and STDERR files to run folder
cd /home/ryank/NGSpipeline_STDOUT_files/
cp $JOB_NAME.o$JOB_ID.$SGE_TASK_ID $JOB_NAME.po$JOB_ID.$SGE_TASK_ID -t $HOME/${date}run${runstoday}_$project



echo "Program has run to last line"
