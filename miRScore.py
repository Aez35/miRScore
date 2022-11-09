#!/usr/local/bin/python3.6

# miRScore is a miRNA scoring tool. 
# Its function is to determine if previously annotated miRNA candidates are of high confidence.
# It does by analyzing the sequences of the candidate miRNA and it's hairpin precursor.
# Fastq files containing reads of the proposed miRNAs must be provided as evidence.
# The first step is to quantify reads of each candidate miRNA and ensure it is present in at least 2 libraries.
# Next, the structure of the candidate miRNA is scored accordingly to criteria in Axtell and Lively et al. 2018.
# The output file contains a list of successful candidates, their scores, and details about each successful candidate.
# Additional outputs include the failed miRNAs and the step in which they failed.

import RNA
import pandas as pd
import RNA
import pysam
from pysam import FastaFile
import argparse
import subprocess
from Bio import SeqIO
import os
import glob
from pathlib import Path
import csv


#Assign arguments for command line
ap = argparse.ArgumentParser()

ap.add_argument("-fd", "--fastqd", required = True, help="Directory containing fastq files")
ap.add_argument("-mat", "--mature", required = True, help="fasta file of mature miRNA sequence")
ap.add_argument("-hp", "--hairpin", required = True, help="fasta file of precursor sequence")

args = vars(ap.parse_args())



### ____________________ Functions ______________________ ###
def get_s_rels(q_rel_start, q_rel_stop, dotbracket):
    # Compute position of miR* given miR location on a dotbracket
 
    # Get a lookup table of paired positions
    pairs = get_pairs(dotbracket)
 
    # determine the 3' end (stop) of miR* from q_rel_start
    s_rel_stop = None
    for attempts, qpos in enumerate(range(q_rel_start, q_rel_stop)):
        if pairs[qpos] is not None:
            s_rel_stop = pairs[qpos] + attempts + 2
            break
   
    # determine the 5' end (start) or miR* from q_rel_stop - 2
    s_rel_start = None
    for attempts, qpos in enumerate(range((q_rel_stop - 2), q_rel_start, -1)):
        if pairs[qpos] is not None:
            s_rel_start = pairs[qpos] - attempts
            break
   
    return s_rel_start, s_rel_stop
   

def get_pairs(dotbracket):
    # Compute positional lookups for a dotbracket structure
    # Note conversion from enumerate's zero-based system to one-based
    pairs = {}
    for i, j in enumerate(dotbracket):
        pairs[(i + 1)] = None
   
    stack = []
    for i, j in enumerate(dotbracket):
        if j == '(':
            stack.append((i + 1))
        elif j == ')':
            bp1 = stack.pop()
            bp2 = i + 1
            pairs[bp1] = bp2
            pairs[bp2] = bp1
    # testing
    #print(dotbracket)
    #print(pairs)
    #sys.exit()
    return pairs
 

def score(mir,seq,mirspos,mirpos,ss):
#Determine score of candidate miRNA based on criteria outlined in Axtell and Meyers et al 2018.
    reason=[]
    mscore = []
    # Test that precursor is less than 300 nt
    if len(seq)>300:
        mscore.append(5)
        reason.append("Precursor > 300 nt")
    else:
        mscore.append(10)
        
    #Test length of miRNA
    if (len(mir)>24) or (len(mir)<20):
        mscore.append(0)
        reason.append("miRNA > 24 or < 20")
    elif (len(mir)>19) and (len(mir)<23):
        mscore.append(10)
    elif (len(mir)>22) and (len(mir)<25):
        mscore.append(5)
        reason.append("23/24 nt miRNA")

    
    #Obtain dotbrackets of mir and mirstar
    mirstar_db=ss[mirspos[0]:mirspos[1]]
    mir_db = ss[mirpos[0]:mirpos[1]]

    #Test that mir/mir* duplex does not contain greater than 5 mismatches
    if (mirstar_db.count(".")>5) or (mir_db.count(".") > 5):
        mscore.append(5)
        reason.append("more than 5 mismatches")
    else:
        mscore.append(10)
        
    if index_of("....", mirstar) == -1 and index_of("...", mir) != -1:
        mscore.append(5)
        reason.append("asymmetric bulge greater than 3")
    elif index_of("....", mirstar) != -1 and index_of("...", mir) == -1:
        mscore.append(5)
        reason.append("asymmetric bulge greater than 3")
    else:
        mscore.append(10)
    return mscore, reason



def index_of(asym, in_list):
    try:
        return in_list.index(asym)
    except ValueError:
        return -1

### _______________________ Code begins __________________________ ###
    
################ Generate BAM Files  ################


print("Building index...")
##Map fastq files to hairpin using bowtie2
subprocess.call(['bowtie-build', "precursors.fa", 'hairpin'])

##Create a list of all fastq files in provided directory.
fastqs=glob.glob(args["fastqd"] + '/*.fastq')
##Make output directory.
subprocess.call(["mkdir","alignments"])

print("Mapping each fastq...")
for file in fastqs:
##map each fastq file in directory to index.
    print("Mapping " + file + "...")
    subprocess.call(['bowtie', "-a", 'hairpin', file , "-S", "alignments/" + Path(file).stem + ".bam"])



################  Check for miRNA Reads  ################
##Initialize dictionary of failed miRNAs and successful candidate miRNAs.
## These dictionaries will be used to generate the final output of the program.
failed = {}
mirnas = {}
mircounts = {}

##Fetch miRNA precursor sequence matching mature candidate miRNA sequence
mir_fa = FastaFile(args["mature"])
hairpin_fa = FastaFile(args["hairpin"])
mir_dict =  SeqIO.index(args["mature"], "fasta")
hp_dict = SeqIO.index(args["hairpin"], "fasta")


#Check that each miRNA matches the precursor it shares a name with.
rev_mir = {}
for seq in mir_dict:
    hp=hp_dict[seq].seq
    mir=mir_dict[seq].seq
    if index_of(mir, hp) == -1:
        failed[seq] = "Hairpin not found" 
    else:
        rev_mir[seq]=mir_dict[seq].seq


##Create list of bam files from bowtie2 output.
bamfiles=glob.glob('alignments/*.bam')
#print(bamfiles)

for bam in bamfiles:
    print("Beginning next bam file...")

    pysam.sort("-o", "alignments/" + Path(bam).stem + "_sorted.bam", bam)
    pysam.index("alignments/" + Path(bam).stem + "_sorted.bam")

    #Extract all sequence reads from bam file and create a list of all reads.
    bamfile = pysam.AlignmentFile("alignments/" + Path(bam).stem + "_sorted.bam", "rb")
    reads=[]
    for read in bamfile.fetch():
        reads.append(read.seq)

    for x in rev_mir:
        reads=[]
        for read in bamfile.fetch(x):
            reads.append(read.seq)
        if reads.count(rev_mir[x]) > 0:
            if mircounts.get(x)== None:
                #If not in {mircounts}, add it.
                #print(x + " not yet in mircounts")
                mircounts[x]= 1
            else:
                mircounts[x] = mircounts[x] + 1

        #else:
            #print("miRNA not found.")
        reads.clear()

for x in mircounts:
    if mircounts[x] > 1:
        mirnas[x] = mir_fa.fetch(x)
    else:
        failed[x] = "present in >2 libraries" 

print(str(len(mirnas)) + " potential miRNAs found.")
print(str(len(failed)) + " miRNAs failed to have reads in 2 or more libraries.")

################  Find miR* and score miRNA  ################

header=["miRNA","miRScore","Precursor Length","mir sequence", "miRNA length","mir* sequence","miR* Length", "Reasons"]
with open('mirscore_output.csv', mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)

    data=[]
    for key in rev_mir:
        #print(key)
        seq = hairpin_fa.fetch(key)
        #print(seq)
        mat = mir_fa.fetch(key)
        #print(mat)
        
        (ss, mfe)=RNA.fold(seq)
        #print(ss)

        mstart = seq.index(mat)+1
        mstop = seq.index(mat)+len(mat)+1
        #print(mstart, mstop)

        mirpos=[mstart,mstop]
        mirspos=get_s_rels(mstart,mstop,ss)
        #print(mirspos)
        #print(seq[mstart:mstop])
        #print(seq[mirspos[0]:mirspos[1]])
        #print(mirspos)
        mirstar=ss[mirspos[0]:mirspos[1]]
        #print(mirstar)
        print("Scoring miRNA "+ key)
        result = score(mat,seq,mirspos,mirpos,ss)
        print(result[1])
        print(result[0])
        if len(result[1])>0:
            save=[key,result[0],len(seq),mat,len(mat),seq[mirspos[0]:mirspos[1]],len(mirstar), result[1]]
        else:
            save=[key,result[0],len(seq),mat,len(mat),seq[mirspos[0]:mirspos[1]],len(mirstar)]
        data.append(save)
    csv_writer.writerow(header)
    for row in data:
        csv_writer.writerow(row)
        
faildata=[]
failheader=["miRNA", "Reason for failing"]
with open('failed_output.csv', mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(failheader)
    for fail in failed:
        failsave=[fail, failed[fail]]
        faildata.append(failsave)
    for rows in faildata:
        csv_writer.writerow(rows)
