#!/bin/bash

## qsub -p -50 MokaSubmission.sh
##Check fastq files are loaded into the GSDdata folder
##Check Mutation reports are loaded into the NextGeneMutationReports folder
##Check that folder with Project number is written as Projectnumber eg Project123 There should be no spaces
##Check number of samples and change task number value accordingly
##Scripts to integrate Moka pipeline on SGE cluster
##Set the working directory to the current one

##Output STDERR and STDOUT to output file place in NGSpipeline_STDOUT_files in Home directory
#$ -j y
#$ -o /home/ryank/NGSpipeline_STDOUT_files
#$ -e /home/ryank/NGSpipeline_STDOUT_files

##Set the shell 
#$ -S /bin/bash

##Set the name of the Job
#$ -N SubmissionJob

##Set the Parallel Environment
#$ -pe threaded 1

#Change directory to where the pre-analysed fastq files are held
cd /home/ryank/GSDdata

##Find the Project number for each fastq file and return a list of Project numbers for these files
list=`find -maxdepth 4 -mindepth 2 -name NGS*_GSD_*fastq | grep -o "Project[0-9]\+" | uniq | grep -o "[0-9]\+"`

##Take the first project 
##firstprojectnumber=${list%[[:space:]][0-9]*}
##Loop through each available project and generate corresponding run folders for each Project
for item in $list:
do 
echo ${item%[':']}
item=${item%[':']}
project="GSD"
project=$project"_"$item

echo $project

cd $HOME
##Make a folder with today's date 

date=$(date +"%d-%m-%y")
echo $date

#Find the number of folders which start with today's date followed by 'run' plus a number and assign it to the variable runstoday
runstoday=`find -mindepth 1 -maxdepth 1 -name ${date}[run]\*[0-9]\* | wc -l`
echo $runstoday
##Make the parent directory
mkdir $HOME/${date}run${runstoday}_$project

##Copy the Moka NGS Pipeline to the new analysis folder
cp /home/ryank/MokaPipeline14-02-2014.sh -t $HOME/${date}run${runstoday}_$project/

##Make child folders for run
mkdir $HOME/${date}run${runstoday}_$project/fastq

mkdir $HOME/${date}run${runstoday}_$project/Alignment_Output

mkdir $HOME/${date}run${runstoday}_$project/VCF_files

mkdir $HOME/${date}run${runstoday}_$project/Annotation_files

##Change directory to where the pre-analysed fastq files are held
cd /home/ryank/GSDdata/Project$item

##find all the fastq files containing the term GSD and mv them to the fastq folder in the run folder for the MokaPipeline 
find -maxdepth 4 -mindepth 3 -name NGS*_GSD_*fastq -exec mv {} $HOME/${date}run${runstoday}_$project/fastq \; 

#Calculate what the upper range for the number of tasks which are required to be submitted
tasknumber=`ls $HOME/${date}run${runstoday}_$project/fastq | wc -l`

tasknumber=$((tasknumber / 2))
echo "tasknumber= $tasknumber"
##Generate a filelist.txt file which can be read by MokaPipeline14-02-2014.sh when it is assigning fastq files to each task
ls -I filelist.txt $HOME/${date}run${runstoday}_$project/fastq/ > $HOME/${date}run${runstoday}_$project/fastq/filelist.txt

# cd to home directory and ssh to apollo node to submit the MokaPipeline14-02-2014.sh job with the associated number of tasks
cd /home/ryank/
ssh ryank@apollo "cd /home/ryank; qsub -p -50 -t 1-$tasknumber -v date=$date,runstoday=$runstoday,project=$project MokaPipeline14-02-2014.sh"

done





