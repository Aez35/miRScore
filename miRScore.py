
#!/usr/local/bin/python3.6

# miRScore is a miRNA scoring tool. 
# Its function is to determine if previously annotated miRNA candidates are of high confidence.
# It does by analyzing the sequences of the candidate miRNA and it's hairpin precursor.
# Fastq files containing reads of the proposed miRNAs must be provided as evidence.
# The first step is to quantify reads of each candidate miRNA and ensure it is present in at least 2 libraries.
# Next, the structure of the candidate miRNA is scored accordingly to criteria in Axtell and Meyers (2018).
# The output file contains a list of successful candidates, their scores, and details about each successful candidate.
# Additional outputs include the failed miRNAs and the step in which they failed.

import RNA
import pandas as pd
import pysam
from pysam import FastaFile
import argparse
import subprocess
from Bio import SeqIO
import os
import glob
from pathlib import Path
import csv
import regex
import sys

#Assign arguments for command line
ap = argparse.ArgumentParser()

ap.add_argument("--fastqd", required = True, help="Directory containing fastq files")
ap.add_argument("--mature", required = True, help="fasta file of mature miRNA sequence")
ap.add_argument("--hairpin", required = True, help="fasta file of precursor sequence")
#ap.add_argument("-mm","--mismatch", help="Allow up to 1 mismatch in miRNA reads")

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
    return pairs
 

def score(mir,seq,mirsspos,mirpos,ss):
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
    mirstar_db=ss[mirsspos[0]:mirsspos[1]]
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
#Index mir/mir* for asymestrical bulges
    try:
        return in_list.index(asym)
    except ValueError:
        return -1

def is_fasta(filename):
#Check if provided file is in fasta format
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

def find_indices_of(char, in_string):
    index = -1
    while True:
        index = in_string.find(char, index + 1)
        if index == -1:
            break
        yield index


### _______________________ Code begins __________________________ ###

#Check if provided hairpin and mature sequences are in fasta format
if is_fasta(args["hairpin"]) == False:
    sys.exit("Error: Hairpin file must be in fasta format.")
if is_fasta(args["mature"]) == False:
    sys.exit("Error: Mature sequence must be in fasta format.")


# Read in candidate miRNA hairpin and precursor
mir_fa = FastaFile(args["mature"])
hairpin_fa = FastaFile(args["hairpin"])


    
#***********  STEP 1: Generate BAM Files  *****************


print("Building index...")
#Map fastq files to hairpin using bowtie2
subprocess.call(['bowtie-build', args["hairpin"], 'hairpin'])

##Create a list of all fastq files in provided directory.
fastqs=glob.glob(args["fastqd"] + '/*.fastq')
if len(fastqs)<2:
    sys.exit("Error: miRScore requires at least two fastq files.")
#Make output directory.
subprocess.call(["mkdir","alignments"])

print("Mapping each fastq...")
for file in fastqs:
#Map each fastq file in directory to index.
    print("Mapping " + file + "...")
    subprocess.call(['bowtie', "-a", 'hairpin', file , "-S", "alignments/" + Path(file).stem + ".bam"])


#***********  STEP 2: Check for miRNA Reads  *****************
#Initialize dictionary of failed miRNAs and successful candidate miRNAs.
## These dictionaries will be used to generate the final output of the program.
failed = {}
mirnas = {}
rev_mir={}
fir_mir={}

#Fetch miRNA precursor sequence matching mature candidate miRNA sequence
mir_dict =  SeqIO.index(args["mature"], "fasta")
hp_dict = SeqIO.index(args["hairpin"], "fasta")

#Check that each miRNA and hairpin contains only A, T, G, C. If not, add to failed or quit.
for x in mir_dict:
    characters = mir_dict[x].seq
    result = all(char in characters for char in 'ATCG')
    if result == False:
        failed[x]="Sequence contained characters besides A, T, G, or C"
    else:
        fir_mir[x]=mir_dict[x].seq


#Check miRNA integrity and mismatches.
for name in fir_mir:
    #Check that miRNA name is found in the hairpin dictionary.
    if name in hp_dict:
        hp=hp_dict[name].seq
        mir=mir_dict[name].seq
        #Check that miRNA does not multimap to hairpin
        if hp.count(mir)>1:
            failed[name] = "miRNA multimaps to hairpin" 
        else:
            #Check for any mismatches in miRNA compared to hairpin.
            m=regex.findall(str(mir) + '{s<=1}', str(hp)) #allow up to 1 error
            if len(m) < 1:
                failed[name] = "Hairpin does not match" 
            else:
                rev_mir[name]=m[0]
    else:
        failed[name] = "Candidate miRNA does not have a hairpin of the same name" 


#Create list of bam files from bowtie2 output.
bamfiles=glob.glob('alignments/*.bam')

#Sort each bamfile and create a dictionary of read counts for each miRNA
mircounts = {}
for bam in bamfiles:
    print("Beginning next bam file...")
    #Sort bam files
    pysam.sort("-o", "alignments/" + Path(bam).stem + "_sorted.bam", bam)
    pysam.index("alignments/" + Path(bam).stem + "_sorted.bam")

    #Read in bamfile in pysam format
    bamfile = pysam.AlignmentFile("alignments/" + Path(bam).stem + "_sorted.bam", "rb")

    #Count the number of reads each miRNA (x) appears in bamfile
    for x in rev_mir:
        reads=[]
        #Get loci of miRNA
        mir=rev_mir[x]
        hp=hp_dict[x].seq

        indices=[]
        for i in find_indices_of(mir, hp):
            indices.append(i)
        if len(indices)==1:
            # Allow 1-2 positional variance in miRNA reads
            mstart=hp.index(mir) - 1
            mstop=mstart+len(mir) + 1
            #Count reads of miRNA in bamfile
            for read in bamfile.fetch(x,mstart,mstop):
                reads.append(read.seq)
            if reads.count(mir) > 0:
                if mircounts.get(x)== None:
                    #If not in {mircounts}, add it.
                    #print(x + " not yet in mircounts")
                    mircounts[x]= 1
                else:
                    mircounts[x] = mircounts[x] + 1
            else:
                if mircounts.get(x) == None:
                    mircounts[x] = 0

            reads.clear()

#If at least one read of the miRNA (x) is found, add it to {mirnas}. Otherwise add it to {failed}.
for x in mircounts:
    if mircounts[x] > 1:
        mirnas[x] = rev_mir[x]
    else:
        failed[x] = "present in <2 libraries" 

print(str(len(mirnas)) + " potential miRNAs found.")
print(str(len(failed)) + " miRNAs failed to have reads in 2 or more libraries.")

#***********  STEP 3: Count miR*, score miR/miR*/hairpin, and write output file  *****************

header=["miRNA","miRScore","Precursor Length","mir sequence", "miRNA length","mir* sequence","miR* Length", "Reasons"]
with open('novel_miRNAs.csv', mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)

    data=[]
    mirstar_counts={}
    for key in mirnas:
        #Hairpin sequence
        seq = hairpin_fa.fetch(key)
        #mature sequence
        mat = rev_mir[key]
        #Vienna RNA Fold 
        (ss, mfe)=RNA.fold(seq)

        #Index mir while allowing 1 nt mismatch
        m=regex.findall(str(mat) + '{s<=1}', str(seq)) # s<=1 allows 1 mismatch

        #Find start and stop of miR
        mstart=seq.index(m[0])
        mstop=mstart+len(mat)

        #Locate start and stop of miR/miR*
        mirpos=[mstart,mstop]
        mirspos=get_s_rels(mstart,mstop,ss)
        mirsspos=[(mirspos[0]-1),(mirspos[1]-1)]
        
        #Count if miR* reads are present in bamfiles.
        reads=[]
        #Get loci of miRNA
        hp=hp_dict[key].seq
        mirstar=hp[mirsspos[0]:mirsspos[1]]
        
        #Get start and stop of miR*
        mirstart=mirsspos[0]
        mirstop=mirsspos[1]

        #Count miR* reads in bamfiles
        bamfiles=glob.glob("alignments/*sorted.bam")
        for bam in bamfiles:
            bamfile = pysam.AlignmentFile(bam, "rb")
            for read in bamfile.fetch(key,mirstart,mirstop):
                reads.append(read.seq)
            if reads.count(mirstar) > 0:
                #print(reads)
                if mirstar_counts.get(key)== None:
                #If not in {mircounts}, add it.
                    print(key + " not yet in mircounts")
                    mirstar_counts[key]= 1
                else:
                    mirstar_counts[key] = mirstar_counts[key] + 1
            else:
                if mirstar_counts.get(key)== None:
                    mirstar_counts[key]= 0
            reads.clear()

        #If counts detected for miR*, score miRNA and save results.
        if mirstar_counts[key]!= None:
            if mirstar_counts[key] > 1:
                    print("Scoring miRNA "+ key)
                    result = score(mat,seq,mirsspos,mirpos,ss)
                    print(result[1])
                    print(result[0])
                    if len(result[1])>0:
                        save=[key,result[0],len(seq),mat,len(mat),seq[mirsspos[0]:mirsspos[1]],len(mirstar), result[1]]
                    else:
                        save=[key,result[0],len(seq),mat,len(mat),seq[mirsspos[0]:mirsspos[1]],len(mirstar)]
                    data.append(save)
            else:
                failed[key]="miR* reads not present"
                

    csv_writer.writerow(header)
    for row in data:
        csv_writer.writerow(row)

#Remove excess bam files and indexes.
args = ('rm', '-rf', 'alignments/*sorted*')
subprocess.call('%s %s %s' % args, shell=True)
args = ('rm', '-rf', '*ebwt')
subprocess.call('%s %s %s' % args, shell=True)

#Write failed output file for any failed miRNAs.
faildata=[]
failheader=["miRNA", "Reason for failing"]
with open('failed_miRNAs.csv', mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(failheader)
    for fail in failed:
        failsave=[fail, failed[fail]]
        faildata.append(failsave)
    for rows in faildata:
        csv_writer.writerow(rows)
#Remove alignments directory.
subprocess.call(["rm","-r","alignments"])
