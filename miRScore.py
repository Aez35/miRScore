#!/usr/bin/env python
#miRScore V0.86

import RNA
import re
import pysam
from pysam import FastaFile
import argparse
import shutil
import subprocess
from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO
import os
import tqdm
import glob
from pathlib import Path
import csv
import sys
import time
import pandas as pd
import numpy as np
startTime = time.time()

#Check required packages
print("Checking required modules...")
req_packs=["bowtie","samtools","RNAfold"]
for pack in req_packs:
        path = shutil.which(pack)
        if path is not None:
            print('Requires package {0} : {1}'.format(pack,path))
        else:
            msg = 'Requires package {0} : Not found!'.format(pack)
            sys.exit(msg)


#Assign arguments for command line
ap = argparse.ArgumentParser()

ap.add_argument('-bamfile',nargs='*',help='One or more BAM alignment files')
ap.add_argument('-fastq',nargs='*',help='One or more fastq alignment files')
ap.add_argument("-mirnas", required = True, help="fasta file of mature miRNA sequence")
ap.add_argument("-precursors", required = True, help="fasta file of precursor sequence")
ap.add_argument("-mm", help="Allow up to 1 mismatch in miRNA reads", nargs='?', const='')
ap.add_argument("-n", help="Results file name")
ap.add_argument("-nostrucvis", help="Do not include StrucVis output",  nargs='?', const='')
ap.add_argument("-threads", help="Specify number of threads for bam2fq")


args = ap.parse_args()

if args.bamfile is not None:
    bams=list(args.bamfile)
    for b in bams:
        if b.lower().endswith(('.bam'))!=True:
            msg = (b + " is not a bamfile. Please revise input.")
            sys.exit(msg)

if args.fastq is not None:
    fastqs=list(args.fastq)
    for f in fastqs:
        if f.lower().endswith(('.fastq','.fq'))!=True:
            msg = (f + " is not a fastq. Please revise input.")
            sys.exit(msg)


#Check paths exist
path="alignments"
isExist = os.path.exists(path)
if isExist==True:
    msg = ("alignment folder exists in directory. Please remove before running again.")
    sys.exit(msg)

path="RNAplots"
isExist = os.path.exists(path)
if isExist==True:
    msg = ("RNAplots folder exists in directory. Please remove before running again.")
    sys.exit(msg)

path="strucVis"
isExist = os.path.exists(path)
if isExist==True:
    msg = ("strucvis folder exists in directory. Please remove before running again.")
    sys.exit(msg)


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
    criteria=[]
    # Test that precursor is less than 300 nt
    if len(seq)>300:
        mscore.append(10)
        reason.append("Precursor > 300 nt")
        criteria.append("N")
    else:
        mscore.append(25)
        criteria.append("Y")
        
    #Test length of miRNA
    if (len(mir)>24) or (len(mir)<20):
        mscore.append(0)
        reason.append("miRNA > 24 or < 20")
        criteria.append("N")
    elif (len(mir)>19) and (len(mir)<23):
        mscore.append(25)
        criteria.append("Y")
    elif (len(mir)>22) and (len(mir)<25):
        mscore.append(25)
        reason.append("23/24 nt miRNA")
        criteria.append("23/24")

    #Obtain dotbrackets of mir and mirstar
    mirstar_db=ss[mirsspos[0]:mirsspos[1]]
    mir_db = ss[mirpos[0]:mirpos[1]]

    #Test that mir/mir* duplex does not contain greater than 5 mismatches
    if (mirstar_db.count(".")>5) or (mir_db.count(".") > 5):
        mscore.append(0)
        reason.append("More than 5 mismatches")
        if (mirstar_db.count("."))>(mir_db.count(".")):
            criteria.append(mirstar_db.count("."))
        else:
            criteria.append(mir_db.count("."))

    else:
        mscore.append(25)
        criteria.append(mir_db.count("."))


    #Test that mir/miR* does not contain an asymmetric bulge greater than 3    
    if index_of("....", mirstar_db) == -1 and index_of("....", mir_db) != -1:
        mscore.append(10)
        reason.append("asymmetric bulge greater than 3")
        criteria.append("1")
    elif index_of("....", mirstar_db) != -1 and index_of("....", mir_db) == -1:
        mscore.append(10)
        reason.append("asymmetric bulge greater than 3")
        criteria.append("1")
    else:
        mscore.append(25)
        criteria.append("0")

    return mscore, reason, criteria

def index_of(asym, in_list):
#Determine asymmetry in mir complex
    try:
        return in_list.index(asym)
    except ValueError:
        return -1

def is_fasta(filename):
#Determine if file is a fasta file
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

def check_fastas(args):
#Check that hairpin and miRNA files are FASTA files
    if is_fasta(args.precursors) == False:
        sys.exit("Error: Hairpin file must be in fasta format.")
    if is_fasta(args.mirnas) == False:
        sys.exit("Error: Mature sequence must be in fasta format.")

def find_indices_of(char, in_string):
#Find index of miRNA in hairpin
    index = -1
    while True:
        index = in_string.find(char, index + 1)
        if index == -1:
            break
        yield index
try:
    from subprocess import DEVNULL  # Python 3.
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

def run(cmd) :
#Run subprocess
    proc = subprocess.call(cmd, shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

def align_fastqs(args, fastqs):
#Build index and map hairpins to fastq files
    print('')
    print("Building index...")
    #Create index of hairpins using bowtie
    subprocess.call(['bowtie-build', 'tmp.prec.fa', 'hairpin'],stdout=DEVNULL,stderr=subprocess.STDOUT)
    #Map reads to hairpin
    print("Mapping each fastq file to hairpin...")
    for file in fastqs:
        if args.mm != None:
            print("Mapping " + file + "with one mismatch...")
            subprocess.call(['bowtie', "-a","-v1","--no-unal", "--norc", '-x','hairpin', file , "-S", "alignments/" + Path(file).stem + ".bam"],stdout=DEVNULL,stderr=subprocess.STDOUT)
        else:
            print("Mapping " + file + "...")
            subprocess.call(['bowtie', "-a","-v0","--no-unal", "--norc", '-x','hairpin', file , "-S", "alignments/" + Path(file).stem + ".bam"],stdout=DEVNULL,stderr=subprocess.STDOUT)


def DNAcheck(dna):
#Check that all characters in provided sequences are ATGCU
    y = dna.upper()
    if re.match("^[ATGCU]+$", y):
        return(2)
    else:
        return(1)
    
def add_values_in_dict(dict, key, list_of_values):
    if key not in dict:
        dict[key] = list()
    dict[key].extend(list_of_values)
    return dict


### _______________________ Code begins __________________________ ###

def main():
    #Create temporary files converting any U's to T's
    #hairpins
    cmd = "cat " + args.precursors + "|tr [:lower:] [:upper:]  | sed -E '/(>)/!s/U/T/g' > tmp.prec.fa"
    run(cmd)
    #mirnas
    cmd = "cat " + args.mirnas + "|tr [:lower:] [:upper:]  | sed -E '/(>)/!s/U/T/g' > tmp.mirs.fa"
    run(cmd)

    #Fetch miRNA precursor sequence matching mature candidate miRNA sequence
    mir_dict =  SeqIO.index("tmp.mirs.fa", "fasta")
    hp_dict = SeqIO.index("tmp.prec.fa", "fasta")

    #Working list of miRNAs that failed and the first criteria in which they failed.
    mir_one={}
    mirna_dict={}
    failed={}

##____________________Checking flags ______________________##
    #Checks length of hairpin and if the miRNA multimaps to it. 
    #If multimapping is detect, miRNA is removed and added to failed dictionary
    #Check that each miRNA contains only A, T, G, C. If not, add to failed or quit.
    print('')
    print("Checking hairpin and miRNA sequences...")
    for mir in mir_dict:
        if mir in hp_dict:
            characters = str(mir_dict[mir].seq)
            result = DNAcheck(characters)
            characters_hp = str(hp_dict[mir].seq)
            result2 = DNAcheck(characters_hp)
            if result == 1 or result2 ==1:
                failed[mir]="Sequence contained characters besides U,A, T, G, or C"
            elif len(characters_hp)<50:
                failed[mir] = "Hairpin is less than 50 basepairs"
            elif characters_hp.count(characters)>1:
                failed[mir] = "miRNA multimaps to hairpin" 
            else:
                mir_one[mir]=mir_dict[mir].seq
        else:
            failed[mir] = "No hairpin sequence detected"


    #If any miRNAs failed for these reasons, quit miRScore.
    if len(failed) != 0:
        for fail in failed:
            print(fail+" failed due to: "+failed[fail])
        sys.exit("Please ammend and rerun miRScore.")


    #Score miRNA criteria for hairpin length, mismatches, and bulges
    for x in tqdm.tqdm(mir_one, desc="Scoring miRNAs",disable=None):
        hp=str(hp_dict[x].seq)
        mir=str(mir_one[x])
        
        if mir in hp:
            mstart=hp.index(mir)
            mstop=(hp.index(mir)+len(mir))
            
            (ss, mfe)=RNA.fold(hp)
            mir_ss=ss[mstart:mstop]
            if mstart == 0:
                sys.exit("Error! " +x + " miRNA start position is the first nucleotide of the hairpin precursor. Cannot determine miRStar as miRNA structure requires a 2 nt 3' overhang. Please ammend or remove and rerun miRScore.")
            else:
                mirspos=get_s_rels(mstart,mstop,ss)

            if mirspos[0] == None:
                failed[x]="Hairpin structure invalid"
            else:
                #Python indexing of mir position
                mirsspos=[(mirspos[0]-1),(mirspos[1]-1)]
                #if start is greater than stop, mirstar fails.
                if mirsspos[0] > mirsspos[1]:
                    flag=["Hairpin structure invalid"]
                    sep = ";"
                    mstart=mstart+1
                    add_values_in_dict(mirna_dict,x,[mir.translate(str.maketrans("tT", "uU")),len(mir_one[x]),mstart,mstop,"NA","NA",(mirsspos[0]+1),mirsspos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),"NA","NA","NA","NA","NA","NA","NA"])
                #If the length of miRstar is greater than 25, mirstar fails.
                elif (mirsspos[1]-mirsspos[0])>(25):
                    flag=["Hairpin structure invalid"]
                    sep = ";"
                    mstart=mstart+1
                    add_values_in_dict(mirna_dict,x,[mir.translate(str.maketrans("tT", "uU")),len(mir_one[x]),mstart,mstop,"NA","NA",(mirsspos[0]+1),mirsspos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),"NA","NA","NA","NA","NA","NA","NA"])
                else:
                    mirstar=hp[mirsspos[0]:mirsspos[1]]
                    mirpos=[mstart,mstop]
                    result = score(mir,hp,mirsspos,mirpos,ss)        
                    mirstar=hp[mirsspos[0]:mirsspos[1]]
                    mstart=mstart+1

                    #Add result from score function to mirna dict
                    sep = ";"
                    if len(result[1]) == 0:
                        flag = ["NA"]
                    else:
                        flag=list(result[1])
                    #If mirna contains both ) and (, it fails for pairing to itself. Else, if the result of scoring is greater than 80, add mirna to dictionary as a pass, and if not fail it. 
                    if mir_ss.count("(") != 0 and mir_ss.count(")") !=0:
                        flag=["miR stucture invalid"]
                        add_values_in_dict(mirna_dict,x,[mir.translate(str.maketrans("tT", "uU")),len(mir_one[x]),mstart,mstop,str(mirstar.translate(str.maketrans("tT", "uU"))),len(mirstar),(mirsspos[0]+1),mirsspos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),result[2][2],result[2][3]])
                    elif sum(result[0]) < 80:
                        add_values_in_dict(mirna_dict,x,[mir.translate(str.maketrans("tT", "uU")),len(mir_one[x]),mstart,mstop,str(mirstar.translate(str.maketrans("tT", "uU"))),len(mirstar),(mirsspos[0]+1),mirsspos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),result[2][2],result[2][3]])
                    else:
                        add_values_in_dict(mirna_dict,x,[mir.translate(str.maketrans("tT", "uU")),len(mir_one[x]),mstart,mstop,str(mirstar.translate(str.maketrans("tT", "uU"))),len(mirstar),(mirsspos[0]+1),mirsspos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Pass",sep.join(flag),result[2][2],result[2][3]])
        else:
            print("Error! "+ x + ": miRNA sequence not found in hairpin sequence. Please ammend and rerun miRScore.")

    #Stop program if all miRNAs failed previous check
    if len(mirna_dict)<1:
        sys.exit("Error! No candidate miRNAs left to score.")


## BEGIN READ CHECK ##

    #If argument -fastq used, map FASTQ to hairpins directly
    if args.fastq is not None:
        #Create a list of all fastq files in provided directory.
        fastqs=list(args.fastq)
        print("Number of Fastq files detected:" + str(len(fastqs)))
        #Make output directory.
        subprocess.call(["mkdir","alignments"])
        align_fastqs(args,fastqs)

    #If argument -bam used, convert BAMs to FASTQs and then map to hairpins
    elif args.bamfile is not None:
        if args.threads is not None:
            subprocess.call(["mkdir","ms_fastqs"])
            print("Converting BAM to FASTQ for mapping to hairpin...")
            bams2convert=list(args.bamfile)
            print("Number of bam files detected:" + str(len(bams2convert)))
            fastqs=[]
            for file in bams2convert:
                cmd = "samtools bam2fq -@ " + args.threads + " " + file + ' > ms_fastqs/' + Path(file).stem + ".fastq"
                run(cmd)
                fastqs.append("ms_fastqs/" + Path(file).stem + ".fastq")
            subprocess.call(["mkdir","alignments"])
            align_fastqs(args,fastqs)
        else:
            subprocess.call(["mkdir","ms_fastqs"])
            print("Converting BAM to FASTQ for mapping to hairpin...")
            bams2convert=list(args.bamfile)
            print("Number of bam files detected:" + str(len(bams2convert)))
            fastqs=[]
            for file in bams2convert:
                cmd = "samtools bam2fq " + file + ' > ms_fastqs/' + Path(file).stem + ".fastq"
                run(cmd)
                fastqs.append("ms_fastqs/" + Path(file).stem + ".fastq")
            subprocess.call(["mkdir","alignments"])
            align_fastqs(args,fastqs)

    bamfiles=glob.glob('alignments/*.bam')

    #Sort bam files
    print("Sorting bamfiles...")

    #Index and sort each bam file
    for bam in bamfiles:
        pysam.sort("-o", "alignments/" + Path(bam).stem + "_s.bam", bam)
        pysam.index("alignments/" + Path(bam).stem + "_s.bam")
        cmd=("rm -rf "+ bam) 
        run(cmd)        

    #Read counting
    mirna_counts={}
    loci_counts={}
    #Merge library files
    cmd="samtools merge -r alignments/merge1.bam alignments/*s.bam"
    run(cmd)
    #Minimum read length of 15
    cmd=("samtools view -e 'qlen>15' -bh alignments/merge1.bam >alignments/merged.bam")
    run(cmd)
    print("Indexing merged bamfile...")
    cmd="samtools index alignments/merged.bam"
    run(cmd)
    cmd="rm -rf alignments/merge1.bam"
    run(cmd)

    #Create a list of read groups in the merged bam file.
    rgs=[]
    bam_input = pysam.AlignmentFile("alignments/merged.bam", 'rb')
    rg=bam_input.header["RG"]
    for group in range(0,len(rg)):
        for values in rg[group].values():
            rgs.append(values)

    #Count number of reads for each miRNA loci
    print("Counting reads for each miRNA locus...")
    for bam in rgs:
        print("Counting reads in " + bam + "...")
        for x in tqdm.tqdm(mirna_dict,disable=None):
            #Set window for miRNA with variance
            if "Hairpin structure invalid" not in mirna_dict[x][11]:
                if mirna_dict[x][2] !=0 and mirna_dict[x][6]!=0:
                    mstart_indx=mirna_dict[x][2]-1
                    mstop_indx=mirna_dict[x][3]+1
                    mirstart_indx=mirna_dict[x][6]-1
                    mirstop_indx=mirna_dict[x][7]+1
                else:
                    print("Error! " + x +" miRNA starts at first position of hairpin precursor. Please ammend and rerun.")
                    sys.exit()

                add_values_in_dict(mirna_counts,x,[str(Path(bam).stem)])

                #MIR COUNTS
                command=("samtools view -r " + Path(bam).stem + " -F 16 -e 'pos>=" +str(mstart_indx)+ " && endpos<=" + str(mstop_indx)+ "' alignments/merged.bam "+ x+" | wc -l")
                p = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
                output = p.stdout.read()
                splt = output.decode()
                mircount=splt.split()
                add_values_in_dict(mirna_counts,x,[mircount[0]])
                
                
                #MIRSTAR counts
                command=("samtools view -r " + Path(bam).stem + " -F 16 -e 'pos>=" +str(mirstart_indx)+ " && endpos<=" + str(mirstop_indx)+ "' alignments/merged.bam " + x+" | wc -l")
                p = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
                output = p.stdout.read()
                splt = output.decode()
                mirscount=splt.split()
                add_values_in_dict(mirna_counts,x,[mirscount[0]])

                #Total counts for miRNA hairpin
                command=("samtools view -r " + Path(bam).stem + " -F 16 alignments/merged.bam "+ x+ " | wc -l")
                p = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
                output = p.stdout.read()
                splt = output.decode()
                totcount=splt.split()
                add_values_in_dict(loci_counts,x,[int(totcount[0])])
                add_values_in_dict(mirna_counts,x,[totcount[0]])

                #Calculate precision
                if int(totcount[0]) != 0:
                    precision="{:.2f}".format((100*((int(mircount[0])+int(mirscount[0]))/int(totcount[0]))))
                    add_values_in_dict(mirna_counts,x,[precision])
                else:
                    precision= 0.00
                    add_values_in_dict(mirna_counts,x,[precision])

    #Create pandas dataframe for writing reads file
    df=pd.DataFrame.from_dict(mirna_counts, orient='index')
    index=[]
    for x in mirna_counts:
        for bam in rgs:
            index.append(x)
    df2=pd.DataFrame(df.values.reshape(-1, 5),columns=['library','mReads','msReads','allReads','precision'])
    df2.index=index
    #Add pass/fail column
    df2['precision'] = df2['precision'].astype('float')
    df2["mReads"] = df2["mReads"].astype(int)
    df2["msReads"] = df2["msReads"].astype(int)
    df2['result'] = np.where((df2['mReads']>0) & (df2['precision']>=75) & (df2['msReads']>0),'Pass','Fail')
    #Save reads file to csv
    if args.n != None:
        df2.to_csv(args.n+'_reads.csv', sep=',', encoding='utf-8', index_label='miRNA')
    else:
        df2.to_csv('reads.csv', sep=',', encoding='utf-8', index_label='miRNA')

    #Determine if miRNA and miRStar are found in multiple libraries. Precision is calculated individually, library by library.
    for x in tqdm.tqdm(mirna_counts,desc="Counting miRNAs",disable=None):
    #Count mir and mirstar reads for each locus
        mircounts=[]
        mirscounts=[]
        mirstar_reads=(mirna_counts[x][2::5])
        mir_reads=(mirna_counts[x][1::5])
        all_reads=(mirna_counts[x][3::5])
        for c in mirstar_reads:
            mirscounts.append(int(c))
        for c in mir_reads:
            mircounts.append(int(c))
        precCount=[]
        p=(mirna_counts[x][4::5])
        for c in p:
            precCount.append(float(c)) 

#Changed from 'and' to 'or' on account that there should be reads detected for both
        if sum(mircounts)==0 or sum(mirscounts)==0:
            mirna_dict[x][10]="Fail"
            if "NA" in mirna_dict[x][11]:
                mirna_dict[x][11]="No miR duplex reads detected"
                if sum(loci_counts[x]) != 0:
                    precision="{:.2f}".format((100*(sum(mircounts)+sum(mirscounts))/sum(loci_counts[x])))
                    add_values_in_dict(mirna_dict,x,[sum(mircounts),sum(mirscounts),sum(loci_counts[x]),len(bamfiles),precision])
                else:
                    precision= 0.00
                    add_values_in_dict(mirna_dict,x,[sum(mircounts),sum(mirscounts),sum(loci_counts[x]),len(bamfiles),precision])
            else:
                    mirna_dict[x][11]+=";No miR duplex reads detected"       
                    if sum(loci_counts[x]) != 0:
                        precision="{:.2f}".format((100*(sum(mircounts)+sum(mirscounts))/sum(loci_counts[x])))
                        add_values_in_dict(mirna_dict,x,[sum(mircounts),sum(mirscounts),sum(loci_counts[x]),len(bamfiles),precision])
                    else:
                        precision= 0.00
                        add_values_in_dict(mirna_dict,x,[sum(mircounts),sum(mirscounts),sum(loci_counts[x]),len(bamfiles),precision]) 
        else:
            #combined miR/miR* reads
            totreads=[]
            #Number of libraries with mir/mir* reads and precision greater than 75
            count=0
            #mirReads
            mr=0
            #mirstarReads
            ms=0
            #Total reads
            all=0
            #Total precision across libraries counted in results
            precSum=0
            #Count libraries that have reads present
            for i in range(0,len(mir_reads)):
                #implement read floor for each library
                tot_c=int(mir_reads[i])+int(mirstar_reads[i])
                totreads.append(tot_c)
                if int(tot_c) > 10 and int(mir_reads[i]) > 0 and int(mirstar_reads[i]) > 0 and precCount[i]>=75:
                    count=count+1
                    precSum=precSum+precCount[i]
                    mr=mr+int(mir_reads[i])
                    ms=ms+int(mirstar_reads[i])
                    all=all+int(all_reads[i])

            #calculate precision of libraries with reads present
            if count >= 1:
                prec=str(round(precSum/count,2))
                add_values_in_dict(mirna_dict,x,[mr,ms,all,count,prec])
            else:
                mirna_dict[x][10]="Fail"
                if sum(loci_counts[x]) != 0:
                    precision="{:.2f}".format((100*(sum(mircounts)+sum(mirscounts))/sum(loci_counts[x])))
                else:
                    precision= 0.00
                #If total count of miR and miR* is less than 10, fail the locus
                 #implement read floor for each library
                floor=10
                #If any miRNAs have a library that meets the read floor, check precision and report it failed for precision.
                if any(num >= floor  for num in totreads):
                    if sum(mircounts)>0 and sum(mirscounts)>0:
                        if "NA" in mirna_dict[x][11]:
                            mirna_dict[x][11]="Precision less than 75%"
                            add_values_in_dict(mirna_dict,x,[sum(mircounts),sum(mirscounts),sum(loci_counts[x]),len(bamfiles),precision])
                        else:
                            mirna_dict[x][11]+=";Precision less than 75%"      
                            add_values_in_dict(mirna_dict,x,[sum(mircounts),sum(mirscounts),sum(loci_counts[x]),len(bamfiles),precision])        
                #Otherwise, if no libraries had at least 10 reads, report it failed for reads.
                else:                    
                    if "NA" in mirna_dict[x][11]:
                        mirna_dict[x][11]="Less than 10 reads"
                        add_values_in_dict(mirna_dict,x,[sum(mircounts),sum(mirscounts),sum(loci_counts[x]),len(bamfiles),precision])
                    else:
                        mirna_dict[x][11]+=";Less than 10 reads"       
                        add_values_in_dict(mirna_dict,x,[sum(mircounts),sum(mirscounts),sum(loci_counts[x]),len(bamfiles),precision])

    #Add failed miRNAs into dictionary
    for mirfail in failed:
        if failed[mirfail]!="Hairpin structure invalid" and failed[mirfail]!="No hairpin sequence detected" and failed[mirfail]!="Sequence contained characters besides U,A, T, G, or C" and failed[mirfail]!="miR not found in hairpin":
            seq = str(hp_dict[mirfail].seq)
                #mature sequence
            mat = str(mir_one[mirfail])
            mstart=seq.index(mat)
            mstop=(seq.index(mat)+len(mat))

            (ss, mfe)=RNA.fold(seq)
            mirspos=get_s_rels(mstart,mstop,ss)
            #Python indexing of mir position
            mirsspos=[(mirspos[0]-1),(mirspos[1]-1)]
            mirstar=hp[mirsspos[0]:mirsspos[1]]
            add_values_in_dict(mirna_dict,mirfail,[mat,len(mat),mstart,mstop,mirstar,len(mirstar),mirsspos[0],mirsspos[1],seq,len(seq),"Fail",str(failed[mirfail]),"NA","NA","NA","NA","NA","NA","NA","NA","NA"])
        else:
                if failed[mirfail]=="miR not found in hairpin":
                        seq= str(hp_dict[mirfail].seq)
                        mat = str(mir_dict[mirfail].seq)
                        add_values_in_dict(mirna_dict,mirfail,[mat,len(mat),"NA","NA","NA","NA","NA","NA",seq,len(seq),"Fail",str(failed[mirfail]),"NA","NA","NA","NA","NA","NA","NA","NA","NA"])
                else:
                        mat = str(mir_dict[mirfail].seq)
                        add_values_in_dict(mirna_dict,mirfail,[mat,len(mat),"NA","NA","NA","NA","NA","NA","NA","NA","Fail",str(failed[mirfail]),"NA","NA","NA","NA","NA","NA","NA","NA","NA"])
    DEVNULL = open(os.devnull, 'w')

    #Write final output file
    if args.n != None:
        from csv import writer
        header=["name","mSeq","mLen","mStart","mStop","msSeq","msLen","msStart",'msStop',"precSeq","precLen","result", "flags","mismatches","bulges",'mReads','msReads',"totReads","numLibraries","precisionScore"]
        with open(args.n+"_miRScore_results.csv", mode="w", newline='') as csv_out:
            csv_writer = writer(csv_out)
            csv_writer.writerow(header)
            for k, v in mirna_dict.items():
                csv_writer.writerow([k, *v])
    else:
        from csv import writer
        header=["name","mSeq","mLen","mStart","mStop","msSeq","msLen","msStart",'msStop',"precSeq","precLen","result", "flags","mismatches","bulges",'mReads','msReads',"totReads","numLibraries","precisionScore"]
        with open("miRScore_results.csv", mode="w", newline='') as csv_out:
            csv_writer = writer(csv_out)
            csv_writer.writerow(header)
            for k, v in mirna_dict.items():
                csv_writer.writerow([k, *v])

    ##############################################
    ##Beginning of mirna alternative analyzer
    print("*************")
    print("Reevaluting failed miRNAs for potential alternatives...")
    print("*************")

    alt_mrnas={}
    pred_dict={}

    for x in mirna_dict:
        if mirna_dict[x][10]=="Fail":
            if mirna_dict[x][11]!="No hairpin sequence detected" and mirna_dict[x][11]!="Hairpin structure invalid":
                    mir=str(mirna_dict[x][0])
                    hp=str(mirna_dict[x][8])
                    cmd= "samtools view -F 16 alignments/merged.bam "+x+"| awk -F '\t' '{print $10}' | sort | uniq -c | sort -nr|awk 'length($2) > 19 && length($2) < 25 {print}'| head -n1 | awk '{print $2'}"
                    p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
                    output = p.stdout.read()
                    splt = output.decode()
                    if splt != '':
                        alt_mrnas[x]=str(splt.strip())
    print("miRNAs predicted from failed miRNAs:")
    print(len(alt_mrnas))

#Score alternative miRNAs
    for x in tqdm.tqdm(alt_mrnas, desc="Scoring alternative miRNAs",disable=None):
        hp=str(hp_dict[x].seq)
        mir=str(alt_mrnas[x])
        
        if mir in hp:
            mstart=hp.index(mir)
            mstop=(hp.index(mir)+len(mir))
            #RNAfold each hairpin
            (ss, mfe)=RNA.fold(hp)
            #Dotbracket structure
            mir_ss=ss[mstart:mstop]

            if mstart == 0:
                mirspos=get_s_rels(mstart+1,mstop+1,ss)
            else:
                mirspos=get_s_rels(mstart,mstop,ss)
            
            mstart=mstart+1
            if mirspos[0] != None:
                #Python indexing of mir position
                mirsspos=[(mirspos[0]-1),(mirspos[1]-1)]
                #if start is greater than stop, mirstar fails.
                if mirsspos[0] < mirsspos[1] and (mirsspos[1]-mirsspos[0])<(25):
                    mirstar=hp[mirsspos[0]:mirsspos[1]]
                    mirpos=[mstart,mstop]
                    result = score(mir,hp,mirsspos,mirpos,ss)        
                    mirstar=hp[mirsspos[0]:mirsspos[1]]
                    #Add result from score function to mirna dict
                    sep = ";"
                    if len(result[1]) == 0:
                        flag = ["NA"]
                    else:
                        flag=list(result[1])
                    #If mirna contains both ) and (, it fails for pairing to itself. Else, if the result of scoring is greater than 80, add mirna to dictionary as a pass, and if not fail it. 
                    if sum(result[0]) > 80:
                        add_values_in_dict(pred_dict,x,[mir.translate(str.maketrans("tT", "uU")),len(alt_mrnas[x]),mstart,mstop,str(mirstar.translate(str.maketrans("tT", "uU"))),len(mirstar),(mirsspos[0]+1),mirsspos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Pass",sep.join(flag),result[2][2],result[2][3]])

    print("miRNAs output from scoring step:")
    print(len(pred_dict))

    pred_counts={}
    predloci_counts={}
    for bam in rgs:
        print("Counting reads in " + bam + "...")
        for x in tqdm.tqdm(pred_dict,disable=None):
            if "Hairpin structure invalid" not in pred_dict[x][11]:
                if pred_dict[x][2] !=0 and pred_dict[x][6]!=0:
                    mstart_indx=pred_dict[x][2]-1
                    mstop_indx=pred_dict[x][3]+1
                    mirstart_indx=pred_dict[x][6]-1
                    mirstop_indx=pred_dict[x][7]+1
                else:
                    mstart_indx=pred_dict[x][2]
                    mstop_indx=pred_dict[x][3]+1
                    mirstart_indx=pred_dict[x][6]
                    mirstop_indx=pred_dict[x][7]+1

                add_values_in_dict(pred_counts,x,[str(Path(bam).stem)])


                #MIR COUNTS
                command=("samtools view -r " + Path(bam).stem + " -F 16 -e 'pos>=" +str(mstart_indx)+ " && endpos<=" + str(mstop_indx)+ "' alignments/merged.bam "+ x+" | wc -l")
                p = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
                output = p.stdout.read()
                splt = output.decode()
                mircount=splt.split()
                add_values_in_dict(pred_counts,x,[mircount[0]])
                
                #MIRSTAR counts
                command=("samtools view -r " + Path(bam).stem + " -F 16 -e 'pos>=" +str(mirstart_indx)+ " && endpos<=" + str(mirstop_indx)+ "' alignments/merged.bam " + x+" | wc -l")
                p = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
                output = p.stdout.read()
                splt = output.decode()
                mirscount=splt.split()
                add_values_in_dict(pred_counts,x,[mirscount[0]])

                #Total counts for miRNA hairpin
                command=("samtools view -r " + Path(bam).stem + " -F 16 alignments/merged.bam "+ x+ " | wc -l")
                p = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
                output = p.stdout.read()
                splt = output.decode()
                totcount=splt.split()
                add_values_in_dict(predloci_counts,x,[int(totcount[0])])
                add_values_in_dict(pred_counts,x,[totcount[0]])

                
                if int(totcount[0]) != 0:
                    precision="{:.2f}".format((100*((int(mircount[0])+int(mirscount[0]))/int(totcount[0]))))
                    add_values_in_dict(pred_counts,x,[precision])
                else:
                    precision= 0.00
                    add_values_in_dict(pred_counts,x,[precision])

    #Create pandas dataframe for writing reads file
    df=pd.DataFrame.from_dict(pred_counts, orient='index')
    index=[]
    for x in pred_counts:
        for bam in bamfiles:
            index.append(x)
    df2=pd.DataFrame(df.values.reshape(-1, 5),columns=['library','mReads','msReads','allReads','precisionScore'])
    df2.index=index
    #Add pass/fail column
    df2['precisionScore'] = df2['precisionScore'].astype('float')
    df2['result'] = np.where(df2['precisionScore']>=75, 'Pass','Fail')
    #Save reads file to csv
    if args.n != None:
        df2.to_csv(args.n+'_alt_reads.csv', sep=',', encoding='utf-8', index_label='miRNA')
    else:
        df2.to_csv('alt_reads.csv', sep=',', encoding='utf-8', index_label='miRNA')

    for x in tqdm.tqdm(list(pred_counts),desc="Counting miRNAs",disable=None):
    #Count mir and mirstar reads for each locus
        mircounts=[]
        mirscounts=[]
        mirstar_reads=(pred_counts[x][2::5])
        mir_reads=(pred_counts[x][1::5])
        all_reads=(pred_counts[x][3::5])
        for c in mirstar_reads:
            mirscounts.append(int(c))
        for c in mir_reads:
            mircounts.append(int(c))
        precCount=[]
        p=(pred_counts[x][4::5])
        for c in p:
            precCount.append(float(c)) 

        if sum(mircounts)==0 or sum(mirscounts)==0:
            del(pred_dict[x])
        else:
            count=0
            mr=0
            ms=0
            all=0
            precSum=0
            #Count count libraries that have reads present
            for i in range(0,len(mir_reads)):
                if int(mir_reads[i]) > 0 and int(mirstar_reads[i]) > 0 and precCount[i]>=75:
                    count=count+1
                    precSum=precSum+precCount[i]
                    mr=mr+int(mir_reads[i])
                    ms=ms+int(mirstar_reads[i])
                    all=all+int(all_reads[i])

            #calculate precision of libraries with reads present
            if count >= 1:
                prec=str(round(precSum/count,2))
                add_values_in_dict(pred_dict,x,[mr,ms,all,count,prec])
            else:
                del(pred_dict[x])

    if args.n != None:
        from csv import writer
        header=["name","mSeq","mLen","mStart","mStop","msSeq","msLen","msStart",'msStop',"precSeq","precLen","result", "flags","mismatches","bulges",'mReads','msReads',"totReads","nlibraries","precision_score"]
        with open(args.n+"_predcted_mirna_results.csv", mode="w", newline='') as csv_out:
            csv_writer = writer(csv_out)
            csv_writer.writerow(header)
            for k, v in pred_dict.items():
                csv_writer.writerow([k, *v])
    else:
        from csv import writer
        header=["name","mSeq","mLen","mStart","mStop","msSeq","msLen","msStart",'msStop',"precSeq","precLen","result", "flags","mismatches","bulges",'mReads','msReads',"totReads","nlibraries","precision_score"]
        with open("alt_mirna_results.csv", mode="w", newline='') as csv_out:
            csv_writer = writer(csv_out)
            csv_writer.writerow(header)
            for k, v in pred_dict.items():
                csv_writer.writerow([k, *v])
                
    #Make directories for post-script images

    subprocess.call(["mkdir","RNAplots"])
    subprocess.call(["mkdir","RNAplots/failed"])
    subprocess.call(["mkdir","RNAplots/passed"])

    if args.nostrucvis == None:
        #Index hairpins for StrucVis and make directories
        subprocess.call(["mkdir","strucVis"])
        subprocess.call(["mkdir","strucVis/failed"])
        subprocess.call(["mkdir","strucVis/passed"])
        cmd='samtools faidx tmp.prec.fa'
        run(cmd)


        for x in tqdm.tqdm(mirna_dict,desc="strucVis",disable=None):
            if mirna_dict[x][11]!="No hairpin sequence detected" and mirna_dict[x][11]!="Hairpin structure invalid":
                #Create strucVis plots for each miRNA
                if "Fail" in mirna_dict[x][10]:
                    hp_len=str(mirna_dict[x][9])
                    cmd="strucVis -b alignments/merged.bam -g tmp.prec.fa -c "+ x + ":1-" + hp_len + " -s plus -p strucVis/failed/" +x+ ".ps -n " + x
                    run(cmd)
                else:
                    hp_len=str(mirna_dict[x][9])
                    cmd="strucVis -b alignments/merged.bam -g tmp.prec.fa -c "+ x + ":1-" + hp_len + " -s plus -p strucVis/passed/" +x+ ".ps -n " + x
                    run(cmd)


    #Create RNAplots depicting miRNA loci
    os.chdir('RNAplots')
    for x in tqdm.tqdm(mirna_dict,desc="RNAplots",disable=None):
        if mirna_dict[x][11]!="Hairpin structure invalid" and mirna_dict[x][11]!="No hairpin sequence detected" and mirna_dict[x][11]!="Sequence contained characters besides U,A, T, G, or C" and mirna_dict[x][11]!="miR not found in hairpin":
                seq = str(mirna_dict[x][8])
                #mature sequence
                mstart=int(mirna_dict[x][2])
                mstop=int(mirna_dict[x][3])
                mirstart=int(mirna_dict[x][6])
                mirstop=int(mirna_dict[x][7])
                #Create ps file from RNAplot
                if "Fail" in mirna_dict[x][10]:
                    os.chdir("failed")
                    f= open(x+".fa", "a")
                    f.write("> " + str(x) + '\n' + seq)
                    f.close()
                    command=("RNAfold " + x+ ".fa | RNAplot --pre '" + str(mstart)+ " " + str(mstop) + " 8 1 0.5 0 omark " + str(mirstart) + " " + str(mirstop) + " 8 0 0.5 1 omark'")
                    subprocess.Popen(command, shell=True)
                    os.chdir('..')                
                else:
                    os.chdir("passed")
                    f= open(x+".fa", "a")
                    f.write("> " + str(x) + '\n' + seq)
                    f.close()
                    command=("RNAfold " + x+ ".fa | RNAplot --pre '" + str(mstart)+ " " + str(mstop) + " 8 1 0.5 0 omark " + str(mirstart) + " " + str(mirstop) + " 8 0 0.5 1 omark'")
                    subprocess.Popen(command, shell=True)
                    os.chdir('..')                
    os.chdir('..')


    #Add legend to passed RNAplots
    files=glob.glob("RNAplots/passed/*.fa")
    for file in files:
        x=Path(file).stem
        mstart=int(mirna_dict[x][2])
        mstop=int(mirna_dict[x][3])
        mirstart=int(mirna_dict[x][6])
        mirstop=int(mirna_dict[x][7])
        cmd="awk -v n=14 -v s='0 0 0 setrgbcolor/Helvetica findfont\\n9 scalefont setfont \\n72 114 moveto \\n(miRNA: "+ x+" ) show' 'NR == n {print s} {print}' RNAplots/passed/"+x+"_ss.ps | awk -v n=14 -v s='/Helvetica findfont\\n9 scalefont setfont \\n72 104 moveto \\n(miRNA locus: "+ str(mstart)+ "-" + str(mstop)+ ") show \\n/Helvetica findfont\\n9 scalefont setfont \\n72 94 moveto \\n(miRNAStar locus: "+ str(mirstart)+ "-" + str(mirstop)+ ") show \\n0 0.5 1 setrgbcolor \\n65 96 4 0 360 arc closepath fill stroke \\n 1 0.5 0 setrgbcolor \\n65 106 4 0 360 arc closepath fill stroke' 'NR==n {print s} {print}'  > RNAplots/passed/"+x+"_plot.ps"
        run(cmd)

    #Add legend to failed RNAplots
    files=glob.glob("RNAplots/failed/*.fa")
    for file in files:
        x=Path(file).stem
        mstart=int(mirna_dict[x][2])
        mstop=int(mirna_dict[x][3])
        mirstart=int(mirna_dict[x][6])
        mirstop=int(mirna_dict[x][7])
        cmd="awk -v n=14 -v s='0 0 0 setrgbcolor/Helvetica findfont\\n9 scalefont setfont \\n72 114 moveto \\n(miRNA: "+ x+" ) show' 'NR == n {print s} {print}' RNAplots/failed/"+x+"_ss.ps | awk -v n=14 -v s='/Helvetica findfont\\n9 scalefont setfont \\n72 104 moveto \\n(miRNA locus: "+ str(mstart)+ "-" + str(mstop)+ ") show \\n/Helvetica findfont\\n9 scalefont setfont \\n72 94 moveto \\n(miRNAStar locus: "+ str(mirstart)+ "-" + str(mirstop)+ ") show \\n0 0.5 1 setrgbcolor \\n65 96 4 0 360 arc closepath fill stroke \\n 1 0.5 0 setrgbcolor \\n65 106 4 0 360 arc closepath fill stroke' 'NR==n {print s} {print}'  > RNAplots/failed/"+x+"_plot.ps"
        run(cmd)

#Remove unecessary files
    cmd="rm -rf RNAplots/passed/*_ss.ps*"
    run(cmd)
    cmd="rm -rf RNAplots/failed/*_ss.ps*"
    run(cmd)
    cmd="rm -rf RNAplots/passed/*.fa"
    run(cmd)
    cmd="rm -rf RNAplots/failed/*.fa"
    run(cmd)
    cmd="rm -rf *.ebwt|rm -rf *.fai|rm -rf tmp* | rm -rf alignmnets/*bai"
    run(cmd)
    cmd="rm -rf ms_fastqs/"
    run(cmd)


    print("Summary")
    print("____________________________________")
    if args.n != None:
        print("Results file:" + args.n + "_miRScore_results")
    print("Number of submitted candidate miRNAs: " + str(len(mir_dict)))
    res_pass = 0
    for key in mirna_dict:
        if mirna_dict[key][10] == "Pass":
            res_pass = res_pass + 1
    print("Total number of miRNAs identified with high confidence: " + str(res_pass) )
    res_fail = 0
    for key in mirna_dict:
        if mirna_dict[key][10] == "Fail":
            res_fail = res_fail + 1
    print("Total number of Failed miRNAs: " + str(res_fail))
    print('')
    executionTime = (time.time() - startTime)
    print('Time to run: ' + str(executionTime))
    print('Run Completed!')


if __name__ == "__main__":
    main()
