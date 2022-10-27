# python3.6 mirscore_step1.py --fastqDir fastqs --mature candidatemiRNAs.fa --hairpin precursors.fa

from tkinter.messagebox import YES
import pandas as pd
import RNA
import seaborn as sns
import pysam
from pysam import FastaFile
import argparse
import subprocess
from Bio import SeqIO
import os
import glob
from pathlib import Path


#Assign arguments for command line input.
ap = argparse.ArgumentParser()
ap.add_argument("-fd", "--fastqDir", help="Directory containing fastq files")
ap.add_argument("-mat", "--mature", help="fasta file of mature miRNA sequence")
ap.add_argument("-hp", "--hairpin", help="fasta file of precursor sequence")
args = vars(ap.parse_args())

# Read in candidate miRNAs and hairpin precursors
mir_fa = FastaFile(args["mature"])
hairpin_fa = FastaFile(args["hairpin"])

#***********  STEP 1.1: Generate BAM File and check for reads  *****************

#Map fastq files to hairpin using bowtie2
subprocess.call(['bowtie2-build', args["hairpin"], 'precursors.index'])

#Create a list of all fastq files in provided directory.
fastqs=glob.glob(args["fastqDir"] + '/*.fastq')
#Make output directory.
subprocess.call(["mkdir","alignments"])

for file in fastqs:
#map each fastq file in directory to index.
    print("Mapping " + file + "...")
    subprocess.call(['bowtie2','-x', 'precursors.index', file , '>' , "alignments/" + Path(file).stem + ".sam"])

#Create list of sam files output from bowtie2.
samfiles=glob.glob('alignments/*.sam')
for sam in samfiles:
    #Convert to bam.
    subprocess.call(["samtools", "view", "-bS", sam],
    stdout = open(os.path.splitext(sam)[0]+".bam", 'w'))

#***********  STEP 1.2: Check for miRNA Reads  *****************
#Initialize dictionary of failed miRNAs and successful candidate miRNAs.
# These dictionaries will be used to generate the final output of the program.
failed = {}
mirnas = {}

#Fetch miRNA precursor sequence matching mature candidate miRNA sequence
mir_dict =  SeqIO.index(args["mature"], "fasta")

#Create list of bam files from bowtie2 output.
bamfiles=glob.glob('alignments/*.bam')

for bam in bamfiles:
    #Extract all sequence reads from bam file and create a list of all reads.
    bamfile = pysam.AlignmentFile("alignments/" + Path(bam).stem + ".bam", "rb")
    reads=[]
    for read in bamfile.fetch():
        reads.append(read.seq)

    for x in mir_dict.values():
    #Determine if each individual  miRNA sequence is detected in bam file. If yes, add to dictionary of mirnas {mirnas}.
    #If not present, add to dictionary of failed mirna candidates {failed}.
        if reads.count(x.seq) > 0:
            #print('reads present')

            if mirnas.get(x.name)== None:
            #If not in {mirnas}, add it.
                print("mirna not yet in dictionary")
                mirnas[x.name]=reads.count(x.seq)
            else:
            #If mirna already exists in {mirnas}, add reads to existing entry.
                print("mirna in dictionary, appending read count")
                mirnas[x.name]=reads.count(x.seq) + mirnas[x.name]
                print(mirnas[x.name])
        else:
        #If no reads detected, add to {failed}.
            failed[x.name]='No reads'
            print("No reads present in " + Path(bam).stem + "... adding to failed")

            if failed.get(x.name) != None and mirnas.get(x.name) != None:
            #Check if miRNA was not detected in at least one of the bam files. If so, remove it from {mirna}.
                print("Found in mirnas, removing now...")
                mirnas.pop(x.name)
            elif failed.get(x.name) != None and mirnas.get(x.name) == None:
                print("Found in failed, but not mirna")
        reads.clear()
