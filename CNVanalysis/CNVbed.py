#!/usr/bin/python
import sys, getopt, os
import pandas as pd
#python CNVbed.py


def bedsplitter():
	bed = pd.read_table('/home/kevin/Documents/NGS_Pipeline/BRCA_FH_Panel/CNVanalysis/BRCA NGS CNV exons plus 40bpJC.txt', header= 0)
	print bed
	Delta = []
	Newstart = []
	Newend = []
	Chromosome = []
	Start = bed.Start
	Stop = bed.End
	print Stop
	
	counter = 0
	for row in Start:
		#Difference = (bed.End[counter] -bed.Start[counter])/50
		#Difference = str(Difference).split('.')
		#Delta.append(Difference[0])
		Newstart.append(bed.Start[counter])
		Chromosome.append(bed.chr[counter])
		while (row + 50 < bed.End[counter]) and (bed.End[counter] -row > 75):
			row = row + 50
			Newend.append(row)
			row += 1
			Newstart.append(row)
			Chromosome.append(bed.chr[counter])
		Newend.append(bed.End[counter]) 
		
		
		counter += 1
		#print Difference
	print len(Newstart)
	print len(Newend)
	print len(Chromosome)
	Chromosome = pd.Series(Chromosome)
	Newstart = pd.Series(Newstart)
	Newend = pd.Series(Newend)
	outputfile = pd.DataFrame(zip(Newstart, Newend),  columns = ["Start", "Stop"], index=[Chromosome])
	outputfile.to_csv(path_or_buf="/home/kevin/Documents/NGS_Pipeline/BRCA_FH_Panel/CNVanalysis/BRCAexonsplitbed", sep='\t')


bedsplitter()
	


