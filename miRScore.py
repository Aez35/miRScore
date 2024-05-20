#!/usr/bin/env python
#miRScore V0.90

import RNA, re, pysam
from pysam import FastaFile
import argparse, shutil, subprocess
from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO
import os, tqdm, glob
from pathlib import Path
import csv, sys, time
import pandas as pd
import numpy as np
from csv import writer
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
ap.add_argument("-mature", required = True, help="fasta file of mature miRNA sequence")
ap.add_argument("-hairpin", required = True, help="fasta file of hairpin precursor sequence")
ap.add_argument("-mm", help="Allow up to 1 mismatch in miRNA reads", nargs='?', const='')
ap.add_argument("-n", help="Results file name")
ap.add_argument("-nostrucvis", help="Do not include StrucVis output",  nargs='?', const='')
ap.add_argument("-threads", help="Specify number of threads for bam2fq")
ap.add_argument("-kingdom", required=True,choices=["plant","animal"],help="Specify animal or plant")
ap.add_argument("-star", help="fasta file of mirna* sequence")
args = ap.parse_args()


print("--------")
print("Options:")
print("     'Mature file' " + args.mature)
if args.star is not None:
    print("     'Star file' " + args.star)
print("     'Hairpin file' " + args.hairpin)
if args.n is not None:
    print("     'Filename' "+ args.n)
if args.mm is not None:
    print("     'Mismatches' Yes")
if args.nostrucvis is not None:
    print("     'No strucvis' Yes")
if args.threads is not None:
    print("     'Threads' " + args.threads)
print("     'Kingdom' "+ args.kingdom)
if args.bamfile is not None:
    print("     'Bamfiles' " +str(args.bamfile))
if args.fastq is not None:
    print("     'Fastqs' " +str(args.fastq))
print("--------")

#Check format of files.
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

DEVNULL = open(os.devnull, 'w')

#______________ FUNCTIONS ________________#

def run(cmd) :
#Run subprocess
    proc = subprocess.call(cmd, shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

def add_values_in_dict(dict, key, list_of_values):
    if key not in dict:
        dict[key] = list()
    dict[key].extend(list_of_values)
    return dict

def process_mirnas(input_file, output_file):
# Execute the shell command using subprocess
    subprocess.run(
        f"cat {input_file} | tr '[:lower:]' '[:upper:]' | sed -E '/(>)/!s/U/T/g' > {output_file}",
        shell=True
    )

def DNAcheck(dna):
#Check that all characters in provided sequences are ATGCU
    y = dna.upper()
    if re.match("^[ATGCU]+$", y):
        return(2)
    else:
        return(1)

def index_of(asym, in_list):
#Determine asymmetry in mir complex
    try:
        return in_list.index(asym)
    except ValueError:
        return -1


def score_animal(mature,star,hp,ss,maturepos,starpos):
#Determine score of candidate matureNA based on criteria outlined in Axtell and Meyers et al 2018.
    reason=[]
    mscore =[]
    criteria=[]
    # Test that precursor is less than 200 nt
    if len(hp)>200:
        mscore.append(10)
        reason.append("Precursor > 200 nt")
        criteria.append("N")
    else:
        mscore.append(20)
        criteria.append("Y")
        
    #Test length of mature
    if (len(mature)>25) or (len(mature)<19):
        mscore.append(0)
        reason.append("Mature miRNA length not met")
        criteria.append("N")
    elif (len(mature)>18) and (len(mature)<26):
        mscore.append(20)
        criteria.append("Y")

    #Test length of star
    if (len(star)>25) or (len(star)<19):
        mscore.append(0)
        reason.append("Star length not met")
        criteria.append("N")
    elif (len(star)>18) and (len(star)<26):
        mscore.append(20)
        criteria.append("Y")

    #Obtain dotbrackets of mature and maturestar
    star_db=ss[starpos[0]:starpos[1]]
    mature_db = ss[maturepos[0]:maturepos[1]]

    #Test that mature/mature* duplex does not contain greater than 5 mismatches
    if (star_db.count(".")>7) or (mature_db.count(".") > 7):
        mscore.append(0)
        reason.append("More than 7 mismatches")
        if (star_db.count("."))>(mature_db.count(".")):
            criteria.append(star_db.count("."))
        else:
            criteria.append(mature_db.count("."))
    elif (star_db.count(".")>5) or (mature_db.count(".") > 5):
        mscore.append(10)
        reason.append("More than 5 mismatches")
        if (star_db.count("."))>(mature_db.count(".")):
            criteria.append(star_db.count("."))
        else:
            criteria.append(mature_db.count("."))
    else:
        mscore.append(20)
        criteria.append(mature_db.count("."))


    #Test that mature/mature* does not contain an asymmetric bulge greater than 3    
    if index_of("....", star_db) == -1 and index_of("....", mature_db) != -1:
        mscore.append(0)
        reason.append("asymmetric bulge greater than 3")
        criteria.append("1")
    elif index_of("....", star_db) != -1 and index_of("....", mature_db) == -1:
        mscore.append(0)
        reason.append("asymmetric bulge greater than 3")
        criteria.append("1")
    else:
        mscore.append(20)
        criteria.append("0")

    return mscore, reason, criteria

def score_plant(mature,star,hp,ss,maturepos,starpos):
#Determine score of candidate matureNA based on criteria outlined in Axtell and Meyers et al 2018.
    reason=[]
    mscore =[]
    criteria=[]
    # Test that precursor is less than 300 nt
    if len(hp)>300:
        mscore.append(10)
        reason.append("Precursor > 300 nt")
        criteria.append("N")
    else:
        mscore.append(20)
        criteria.append("Y")
        
    #Test length of miRNA
    if (len(mature)>24) or (len(mature)<20):
        mscore.append(0)
        reason.append("Mature miRNA length not met")
        criteria.append("N")
    elif (len(mature)>19) and (len(mature)<23):
        mscore.append(20)
        criteria.append("Y")
    elif (len(mature)>22) and (len(mature)<25):
        mscore.append(20)
        reason.append("23/24 nt miRNA")
        criteria.append("23/24")

    #Test length of miRNA
    if (len(star)>24) or (len(star)<20):
        mscore.append(0)
        reason.append("Star length not met")
        criteria.append("N")
    elif (len(star)>19) and (len(star)<23):
        mscore.append(20)
        criteria.append("Y")
    elif (len(star)>22) and (len(star)<25):
        mscore.append(20)
        reason.append("23/24 nt miRNA star")
        criteria.append("23/24")

    #Obtain dotbrackets of mature and maturestar
    star_db=ss[starpos[0]:starpos[1]]
    mature_db = ss[maturepos[0]:maturepos[1]]

    #Test that mature/mature* duplex does not contain greater than 5 mismatches
    if (star_db.count(".")>5) or (mature_db.count(".") > 5):
        mscore.append(0)
        reason.append("More than 5 mismatches")
        if (star_db.count("."))>(mature_db.count(".")):
            criteria.append(star_db.count("."))
        else:
            criteria.append(mature_db.count("."))
    else:
        mscore.append(20)
        criteria.append(mature_db.count("."))


    #Test that mature/mature* does not contain an asymmetric bulge greater than 3    
    if index_of("....", star_db) == -1 and index_of("....", mature_db) != -1:
        mscore.append(0)
        reason.append("asymmetric bulge greater than 3")
        criteria.append("1")
    elif index_of("....", star_db) != -1 and index_of("....", mature_db) == -1:
        mscore.append(0)
        reason.append("asymmetric bulge greater than 3")
        criteria.append("1")
    else:
        mscore.append(20)
        criteria.append("0")

    return mscore, reason, criteria

def get_s_rels(q_rel_start, q_rel_stop, dotbracket):
    """
    Determine the start and stop of miR* for a given miRNA and hairpin precursor
    """
    pairs = get_pairs(dotbracket)
 
    # determine the 3' end of miR* from q_rel_start
    s_rel_stop = None
    for attempts, qpos in enumerate(range(q_rel_start, q_rel_stop)):
        if pairs[qpos] is not None:
            s_rel_stop = pairs[qpos] + attempts + 2
            break
   
    # determine the 5' end or miR* from q_rel_stop - 2
    s_rel_start = None
    for attempts, qpos in enumerate(range((q_rel_stop - 2), q_rel_start, -1)):
        if pairs[qpos] is not None:
            s_rel_start = pairs[qpos] - attempts
            break
    return s_rel_start, s_rel_stop
   
def get_pairs(dotbracket):
    """
    Compute positional lookups for a dotbracket structure
    """
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

def align_fastqs(args, fastqs):
#Build index and map hairpins to fastq files
    print('')
    print("Building index...")
    #Create index of hairpins using bowtie
    subprocess.call(['bowtie-build', 'tmp.hairpin.fa', 'hairpin'],stdout=DEVNULL,stderr=subprocess.STDOUT)
    #Map reads to hairpin
    print("Mapping each fastq file to hairpin...")
    for file in fastqs:
        if args.mm != None:
            print("Mapping " + file + "with one mismatch...")
            subprocess.call(['bowtie', "-a","-v1","--no-unal", "--norc", '-x','hairpin', file , "-S", "alignments/" + Path(file).stem + ".bam"],stdout=DEVNULL,stderr=subprocess.STDOUT)
        else:
            print("Mapping " + file + "...")
            subprocess.call(['bowtie', "-a","-v0","--no-unal", "--norc", '-x','hairpin', file , "-S", "alignments/" + Path(file).stem + ".bam"],stdout=DEVNULL,stderr=subprocess.STDOUT)

def get_counts(bam, x, start_indx=None, stop_indx=None):
    if start_indx is None or stop_indx is None:
        command = f"samtools view -r {Path(bam).stem} -F 16 alignments/merged.bam {x} | wc -l"
    else:
        command = f"samtools view -r {Path(bam).stem} -F 16 -e 'pos>={start_indx} && endpos<={stop_indx}' alignments/merged.bam {x} | wc -l"

    output = subprocess.check_output(command, shell=True, text=True).strip()
    count = int(output)
    return count

def predict_alt_mirnas(mirna_dict):
    alt_mirnas = {}
    for x in mirna_dict:
        if mirna_dict[x][10] == "Fail":
            if mirna_dict[x][11] != "No hairpin sequence detected" and mirna_dict[x][11] != "Hairpin structure invalid":
                mir = str(mirna_dict[x][0])
                hp = str(mirna_dict[x][8])
                cmd = f"samtools view -F 16 alignments/merged.bam {x} | awk -F '\t' '{{print $10}}' | sort | uniq -c | sort -nr | awk 'length($2) > 19 && length($2) < 25 {{print}}' | head -n1 | awk '{{print $2}}'"
                p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
                output = p.stdout.read()
                splt = output.decode()
                if splt != '':
                    alt_mirnas[x] = str(splt.strip())
    print("miRNAs predicted from failed miRNAs:")
    print(len(alt_mirnas))
    return alt_mirnas

def score_alternative_mirnas(alt_mrnas, hp_dict):
    pred_dict = {}

    for x in tqdm.tqdm(alt_mrnas, desc="Scoring alternative miRNAs", disable=None):
        hp = str(hp_dict[x].seq)
        mature = str(alt_mrnas[x])

        if mature in hp:
            mstart = hp.index(mature)
            mstop = (hp.index(mature) + len(mature))
            # RNAfold each hairpin
            (ss, mfe) = RNA.fold(hp)
            # Dotbracket structure

            if mstart == 0:
                mirspos = get_s_rels(mstart + 1, mstop + 1, ss)
            else:
                mirspos = get_s_rels(mstart, mstop, ss)

            mstart = mstart + 1
            if mirspos[0] is not None:
                # Python indexing of mir position
                mirsspos = [(mirspos[0] - 1), (mirspos[1])]
                # if start is greater than stop, mirstar fails.
                if mirsspos[0] < mirsspos[1] and (mirsspos[1] - mirsspos[0]) < (25):
                    mirstar = hp[mirsspos[0]:mirsspos[1]]
                    mirpos = [mstart, mstop]
                    if args.kingdom == "animal":
                        result = score_animal(mature, mirstar, hp,ss, mirpos, mirsspos)
                    if args.kingdom == "plant":
                        result = score_plant(mature, mirstar, hp,ss, mirpos, mirsspos)
                    # Add result from score function to mirna dict
                    sep = ";"
                    if len(result[1]) == 0:
                        flag = ["NA"]
                    else:
                        flag = list(result[1])
                    # If mirna contains both ) and (, it fails for pairing to itself. Else, if the result of scoring is greater than 80, add mirna to dictionary as a pass, and if not fail it.
                    if sum(result[0]) > 81:
                        add_values_in_dict(pred_dict, x, [mature.translate(str.maketrans("tT", "uU")),
                                                          len(alt_mrnas[x]), mstart, mstop,
                                                          str(mirstar.translate(str.maketrans("tT", "uU"))),
                                                          len(mirstar), (mirsspos[0] + 1), mirsspos[1],
                                                          hp.translate(str.maketrans("tT", "uU")), len(hp),
                                                          "Pass", sep.join(flag), result[2][3], result[2][4]])
    return pred_dict
### _______________________ Code begins __________________________ ###
def main():
    #Create temporary files converting any U's to T's
    #hairpins
    process_mirnas(args.hairpin, "tmp.hairpin.fa")
    #miRNA
    process_mirnas(args.mature, "tmp.mature.fa")
    #miRNA*
    if args.star is not None:
        process_mirnas(args.star, "tmp.star.fa")

    #Fetch miRNA precursor sequence matching mature candidate miRNA sequence
    mature_dict =  SeqIO.index("tmp.mature.fa", "fasta")
    hp_dict = SeqIO.index("tmp.hairpin.fa", "fasta")
    if args.star is not None:
        star_dict =  SeqIO.index("tmp.star.fa", "fasta")

    print("miRNAs submitted: " + str(len(mature_dict.values())))

    #Working list of miRNAs that failed and the first criteria in which they failed.
    initial={}
    mirna_dict={}
    failed={}
    #Checks length of hairpin and if the miRNA multimaps to it. 
    #If multimapping is detect, miRNA is removed and added to failed dictionary
    #Check that each miRNA contains only A, T, G, C. If not, add to failed or quit.
    y=[]
    print('')
    print("Checking hairpin and miRNA sequences...")
    for mir in mature_dict:
        mature=str(mature_dict[mir].seq)
        if args.star != None:
            #print("Finding star in dictionary")
            if mir in star_dict:
                star=str(star_dict[mir].seq)
                #print(star)
            else:
                sys.exit( "Error! " + mir + " not found in star file.")
        if mir in hp_dict:
            hp=str(hp_dict[mir].seq)
        else:
            sys.exit( "Error! " + mir + " not found in hairpin file.")
        
        #Check mature, star, and hp sequences for ATGCU
        hp_check=DNAcheck(hp)
        if args.star != None:
            star_check=DNAcheck(star)
        else:
            star_check=2
        mature_check=DNAcheck(mature)
        if hp_check == 1 or star_check ==1 or mature_check == 1:
            sys.exit( "Error!" + mir + " has characters besides A, T, G, C, or U in the sequence. Please ammend or remove this sequence and rerun.")
            
        if args.star is not None:
            #Check hairpin length and if mature or star multimaps to hairpin
            if mature not in hp:
                failed[mir] = "miRNA not found in hairpin sequence"
            if star not in hp:
                failed=[mir]='miRNA not found in hairpin sequence'
            elif len(hp) < 50:
                failed[mir] = "Hairpin is less than 50 basepairs" 
            elif hp.count(mature) > 1:
                failed[mir] = "miRNA multimaps to hairpin" 
            elif hp.count(star)>1:
                failed[mir] = "miRNA multimaps to hairpin" 
            else:
                initial[mir]=str(mature_dict[mir].seq)
        else:
            #Check hairpin length and if mature or star multimaps to hairpin
            if mature not in hp:
                failed[mir] = "miRNA not found in hairpin sequence"
            elif len(hp) < 50:
                failed[mir] = "Hairpin is less than 50 basepairs" 
            elif hp.count(mature) > 1:
                failed[mir] = "miRNA multimaps to hairpin" 
            else:
                initial[mir]=str(mature_dict[mir].seq)        
    print('miRNAs failed: ' + str(len(failed.values())))

    for x in tqdm.tqdm(initial, desc="Scoring miRNAs",disable=None):
        mature=str(mature_dict[x].seq)
        if args.star != None:
            star=str(star_dict[x].seq)
        hp=str(hp_dict[x].seq)

        maturestarti=hp.index(mature)
        maturestop=(maturestarti+len(mature))
        maturestart=maturestarti+1
        if maturestarti == 0:
            sys.exit("Error! " +x + " miRNA start position is the first nucleotide of the hairpin precursor. Cannot determine star as miRNA structure requires a 2 nt 3' overhang. Please ammend or remove and rerun miRScore.")
        (ss, mfe)=RNA.fold(hp)
        mir_ss=ss[maturestarti:maturestop]
        #Determine star position and sequence.
        if args.star is not None:
            starstart=hp.index(star)
            starstop=(starstart+len(star))
            starpos=get_s_rels(maturestarti,maturestop,ss)
            starpos=[(starpos[0]-1),(starpos[1]-1)]          
            if starstart!=starpos[0] or starstop!=starpos[1]:
                if abs(starstart-starpos[0])>1 or abs(starstop-starpos[1])>1:
                    flag=["No 2nt 3' overhang"]
                    print(flag)
                    sep = ";"
                    add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(str(mature)),maturestart,maturestop,str(star),len(str(star)),starstart,starstop,hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),"NA","NA","NA","NA","NA","NA","NA"])
                elif args.kingdom == "animal":
                    result=score_animal(mature,star,hp,ss,maturepos,starpos)
                elif args.kingdom == "plant": 
                    result=score_plant(mature,star,hp,ss,maturepos,starpos)
                #Add result from score function to mirna dict
                sep = ";"
                if len(result[1]) == 0:
                    flag = ["NA"]
                else:
                    flag=list(result[1])
                    #If mirna contains both ) and (, it fails for pairing to itself. Else, if the result of scoring is greater than 80, add mirna to dictionary as a pass, and if not fail it. 
                if mir_ss.count("(") != 0 and mir_ss.count(")") !=0:
                    flag=["mature stucture invalid"]
                    add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,str(star.translate(str.maketrans("tT", "uU"))),len(star),(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),result[2][3],result[2][4]])
                elif sum(result[0]) < 81:
                    add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,str(star.translate(str.maketrans("tT", "uU"))),len(star),(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),result[2][3],result[2][4]])
                else:
                        add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,str(star.translate(str.maketrans("tT", "uU"))),len(star),(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Pass",sep.join(flag),result[2][3],result[2][4]])
            else:
                star=hp[starpos[0]:starpos[1]]
                maturepos=[maturestarti,maturestop]
                if args.kingdom == "animal":
                    result=score_animal(mature,star,hp,ss,maturepos,starpos)
                elif args.kingdom == "plant": 
                    result=score_plant(mature,star,hp,ss,maturepos,starpos)
                #Add result from score function to mirna dict
                sep = ";"
                if len(result[1]) == 0:
                    flag = ["NA"]
                else:
                    flag=list(result[1])
                    #If mirna contains both ) and (, it fails for pairing to itself. Else, if the result of scoring is greater than 80, add mirna to dictionary as a pass, and if not fail it. 
                if mir_ss.count("(") != 0 and mir_ss.count(")") !=0:
                    flag=["mature stucture invalid"]
                    add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,str(star.translate(str.maketrans("tT", "uU"))),len(star),(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),result[2][3],result[2][4]])
                elif sum(result[0]) < 81:
                    add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,str(star.translate(str.maketrans("tT", "uU"))),len(star),(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),result[2][3],result[2][4]])
                else:
                        add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,str(star.translate(str.maketrans("tT", "uU"))),len(star),(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Pass",sep.join(flag),result[2][3],result[2][4]])
        else:
            starpos=get_s_rels(maturestarti,maturestop,ss)
            if starpos[0] == None:
                failed[x]="Hairpin structure invalid"
            else:
                #Python indexing of mir position
                starpos=[(starpos[0]-1),(starpos[1]-1)]
                if starpos[0] > starpos[1]:
                    flag=["Hairpin structure invalid"]
                    sep = ";"
                    add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,"NA","NA",(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),"NA","NA","NA","NA","NA","NA","NA"])
                else:
                    star=hp[starpos[0]:starpos[1]]
                    maturepos=[maturestarti,maturestop]
                    if args.kingdom == "animal":
                        result=score_animal(mature,star,hp,ss,maturepos,starpos)
                    elif args.kingdom == "plant": 
                        result=score_plant(mature,star,hp,ss,maturepos,starpos)
                    #Add result from score function to mirna dict
                    sep = ";"
                    if len(result[1]) == 0:
                        flag = ["NA"]
                    else:
                        flag=list(result[1])
                        #If mirna contains both ) and (, it fails for pairing to itself. Else, if the result of scoring is greater than 80, add mirna to dictionary as a pass, and if not fail it. 
                    if mir_ss.count("(") != 0 and mir_ss.count(")") !=0:
                        flag=["mature stucture invalid"]
                        y.append(x)
                        add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,str(star.translate(str.maketrans("tT", "uU"))),len(star),(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),result[2][3],result[2][4]])
                    elif sum(result[0]) < 81:
                        y.append(x)
                        add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,str(star.translate(str.maketrans("tT", "uU"))),len(star),(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Fail",sep.join(flag),result[2][3],result[2][4]])
                    else:
                        add_values_in_dict(mirna_dict,x,[mature.translate(str.maketrans("tT", "uU")),len(initial[x]),maturestart,maturestop,str(star.translate(str.maketrans("tT", "uU"))),len(star),(starpos[0]+1),starpos[1],hp.translate(str.maketrans("tT", "uU")),len(hp),"Pass",sep.join(flag),result[2][3],result[2][4]])
    #Stop program if all miRNAs failed previous check
    if len(mirna_dict)<1:
        sys.exit("Error! No candidate miRNAs left to score.")
    print("miRNAs that failed scoring step: " + str(len(y)))

    ## BEGIN READ CHECK ##
    DEVNULL = open(os.devnull, 'w')

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
        threads_option = ["-@", args.threads] if args.threads is not None else []
        subprocess.call(["mkdir", "ms_fastqs"])
        print("Converting BAM to FASTQ for mapping to hairpin...")
        bams2convert = list(args.bamfile)
        print("Number of bam files detected:" + str(len(bams2convert)))
        fastqs = []
        for file in bams2convert:
            cmd = ["samtools", "bam2fq"] + threads_option + [file, '>', f'ms_fastqs/{Path(file).stem}.fastq']
            run(" ".join(cmd))
            fastqs.append(f"ms_fastqs/{Path(file).stem}.fastq")
        subprocess.call(["mkdir", "alignments"])
        align_fastqs(args, fastqs)
    
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
    # Merge library files, filter by minimum read length, index merged BAM file, and remove temporary file
    run("samtools merge -r alignments/merge1.bam alignments/*s.bam")
    run("samtools view -e 'qlen>15' -bh alignments/merge1.bam > alignments/merged.bam")
    print("Indexing merged BAM file...")
    run("samtools index alignments/merged.bam")
    run("rm -rf alignments/merge1.bam")

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

                mir_count = get_counts(bam, x, mstart_indx, mstop_indx)
                mirs_count = get_counts(bam, x, mirstart_indx, mirstop_indx)
                tot_count = get_counts(bam, x)

                add_values_in_dict(mirna_counts, x, [mir_count, mirs_count,tot_count])
                
                if int(tot_count) != 0:
                    precision="{:.2f}".format((100*((int(mir_count)+int(mirs_count))/int(tot_count))))
                    add_values_in_dict(mirna_counts,x,[precision])
                else:
                    precision= 0.00
                    add_values_in_dict(mirna_counts,x,[precision])

    #Create pandas dataframe for writing reads file
    df = pd.DataFrame.from_dict(mirna_counts, orient='index')
    index = [x for x in mirna_counts for _ in bamfiles]

    # Reshape DataFrame and set column names
    df2 = pd.DataFrame(df.values.reshape(-1, 5), index=index,
                    columns=['library', 'mReads', 'msReads', 'allReads', 'precision'])

    #Add pass/fail column
    df2['precision'] = df2['precision'].astype('float')
    df2["mReads"] = df2["mReads"].astype(int)
    df2["msReads"] = df2["msReads"].astype(int)
    df2['result'] = df2.apply(lambda row: 'Pass' if row['mReads'] > 0 and row['precision'] >= 75 and row['msReads'] > 0 else 'Fail', axis=1)

    # Save DataFrame to CSV
    filename = args.n + '_reads.csv' if args.n is not None else 'reads.csv'
    df2.to_csv(filename, sep=',', encoding='utf-8', index_label='miRNA')

    #Determine if miRNA and miRStar are found in multiple libraries. Precision is calculated individually, library by library.
    for x in tqdm.tqdm(mirna_counts, desc="Counting miRNAs", disable=None):
        mirstar_reads = list(map(int, mirna_counts[x][2::5]))
        mir_reads = list(map(int, mirna_counts[x][1::5]))
        all_reads = list(map(int, mirna_counts[x][3::5]))
        precCount = list(map(float, mirna_counts[x][4::5]))

        #Changed from 'and' to 'or' on account that there should be reads detected for both
        if sum(mir_reads)==0 or sum(mirstar_reads)==0:
            mirna_dict[x][10]="Fail"
            if "NA" in mirna_dict[x][11]:
                mirna_dict[x][11]="No mature or star reads detected"
                if sum(all_reads) != 0:
                    precision="{:.2f}".format((100*(sum(mir_reads)+sum(mirstar_reads))/sum(all_reads)))
                    add_values_in_dict(mirna_dict,x,[sum(mir_reads),sum(mirstar_reads),sum(all_reads),len(bamfiles),precision])
                else:
                    precision= 0.00
                    add_values_in_dict(mirna_dict,x,[sum(mir_reads),sum(mirstar_reads),sum(all_reads),len(bamfiles),precision])
            else:
                    mirna_dict[x][11]+=";No mature or star reads detected"       
                    if sum(all_reads) != 0:
                        precision="{:.2f}".format((100*(sum(mir_reads)+sum(mirstar_reads))/sum(all_reads)))
                        add_values_in_dict(mirna_dict,x,[sum(mir_reads),sum(mirstar_reads),sum(all_reads),len(bamfiles),precision])
                    else:
                        precision= 0.00
                        add_values_in_dict(mirna_dict,x,[sum(mir_reads),sum(mirstar_reads),sum(all_reads),len(bamfiles),precision]) 
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
                if sum(all_reads) != 0:
                    precision="{:.2f}".format((100*(sum(mir_reads)+sum(mirstar_reads))/sum(all_reads)))
                else:
                    precision= 0.00
                #If total count of miR and miR* is less than 10, fail the locus
                #implement read floor for each library
                floor=10
                #If any miRNAs have a library that meets the read floor, check precision and report it failed for precision.
                if any(num >= floor  for num in totreads):
                    if sum(mir_reads)>0 and sum(mirstar_reads)>0:
                        if "NA" in mirna_dict[x][11]:
                            mirna_dict[x][11]="Precision less than 75%"
                            add_values_in_dict(mirna_dict,x,[sum(mir_reads),sum(mirstar_reads),sum(all_reads),len(bamfiles),precision])
                        else:
                            mirna_dict[x][11]+=";Precision less than 75%"      
                            add_values_in_dict(mirna_dict,x,[sum(mir_reads),sum(mirstar_reads),sum(all_reads),len(bamfiles),precision])        
                #Otherwise, if no libraries had at least 10 reads, report it failed for reads.
                else:                    
                    if "NA" in mirna_dict[x][11]:
                        mirna_dict[x][11]="Less than 10 reads"
                        add_values_in_dict(mirna_dict,x,[sum(mir_reads),sum(mirstar_reads),sum(all_reads),len(bamfiles),precision])
                    else:
                        mirna_dict[x][11]+=";Less than 10 reads"       
                        add_values_in_dict(mirna_dict,x,[sum(mir_reads),sum(mirstar_reads),sum(all_reads),len(bamfiles),precision])

    fails=['Hairpin is less than 50 basepairs','miRNA multimaps to hairpin','miRNA not found in hairpin sequence','Sequence contained characters besides U,A, T, G, or C',"No hairpin sequence detected","Hairpin structure invalid"]
    #Add failed miRNAs into dictionary
    for mirfail in failed:
        if failed[mirfail] not in fails:
            seq = str(hp_dict[mirfail].seq)
                #mature sequence
            mat = str(mature_dict[mirfail])
            mstart=seq.index(mat)
            mstop=(seq.index(mat)+len(mat))

            (ss, mfe)=RNA.fold(seq)
            mirspos=get_s_rels(mstart,mstop,ss)
            #Python indexing of mir position
            mirsspos=[(mirspos[0]-1),(mirspos[1]-1)]
            mirstar=hp[mirsspos[0]:mirsspos[1]]
            add_values_in_dict(mirna_dict,mirfail,[mat,len(mat),mstart,mstop,mirstar,len(mirstar),mirsspos[0],mirsspos[1],seq,len(seq),"Fail",str(failed[mirfail]),"NA","NA","NA","NA","NA","NA","NA","NA","NA"])
        else:
                if failed[mirfail]=="mature miRNA not found in hairpin sequence":
                        seq= str(hp_dict[mirfail].seq)
                        mat = str(mature_dict[mirfail].seq)
                        add_values_in_dict(mirna_dict,mirfail,[mat,len(mat),"NA","NA","NA","NA","NA","NA",seq,len(seq),"Fail",str(failed[mirfail]),"NA","NA","NA","NA","NA","NA","NA","NA","NA"])
                else:
                        mat = str(mature_dict[mirfail].seq)
                        add_values_in_dict(mirna_dict,mirfail,[mat,len(mat),"NA","NA","NA","NA","NA","NA","NA","NA","Fail",str(failed[mirfail]),"NA","NA","NA","NA","NA","NA","NA","NA","NA"])
    DEVNULL = open(os.devnull, 'w')

    #Write final output file
    header = ["name", "mSeq", "mLen", "mStart", "mStop", "msSeq", "msLen", "msStart", "msStop", "precSeq", "precLen", "result", "flags", "mismatches", "bulges", "mReads", "msReads", "totReads", "numLibraries", "precision"]
    output_file = args.n + "_miRScore_results.csv" if args.n else "miRScore_results.csv"

    with open(output_file, mode="w", newline='') as csv_out:
        csv_writer = writer(csv_out)
        csv_writer.writerow(header)
        for k, v in mirna_dict.items():
            csv_writer.writerow([k, *v])

    ##############################################
    ##Beginning of mirna alternative analyzer
    print("*************")
    print("Reevaluting failed miRNAs for potential alternatives...")
    print("*************")

    alt_mirnas=predict_alt_mirnas(mirna_dict)
    pred_dict=score_alternative_mirnas(alt_mirnas,hp_dict)

    print("miRNAs predicted from failed miRNA loci:")
    print(len(alt_mirnas))

    pred_counts={}
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

                mir_count = get_counts(bam, x, mstart_indx, mstop_indx)
                mirs_count = get_counts(bam, x, mirstart_indx, mirstop_indx)
                tot_count = get_counts(bam, x)
                    
                add_values_in_dict(pred_counts, x, [mir_count, mirs_count,tot_count])
                
                if int(tot_count) != 0:
                    precision="{:.2f}".format((100*((int(mir_count)+int(mirs_count))/int(tot_count))))
                    add_values_in_dict(pred_counts,x,[precision])
                else:
                    precision= 0.00
                    add_values_in_dict(pred_counts,x,[precision])

    # Create DataFrame from dictionary
    df = pd.DataFrame.from_dict(pred_counts, orient='index')

    # Repeat the index for each BAM file
    index = [x for x in pred_counts for _ in bamfiles]

    # Reshape DataFrame and set column names
    df2 = pd.DataFrame(df.values.reshape(-1, 5), index=index,
                    columns=['library', 'mReads', 'msReads', 'allReads', 'precisionScore'])

    # Convert 'precisionScore' to float and add 'result' column
    df2['precisionScore'] = df2['precisionScore'].astype(float)
    df2['result'] = df2['precisionScore'].apply(lambda x: 'Pass' if x >= 75 else 'Fail')

    # Save DataFrame to CSV
    filename = args.n + '_alt_reads.csv' if args.n is not None else 'alt_reads.csv'
    df2.to_csv(filename, sep=',', encoding='utf-8', index_label='miRNA')

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

    output_file = args.n + "_alt_mirna_results.csv" if args.n else "alt_mirna_results.csv"

    header = ["name", "mSeq", "mLen", "mStart", "mStop", "msSeq", "msLen", "msStart", 'msStop',
            "precSeq", "precLen", "result", "flags", "mismatches", "bulges", 'mReads', 'msReads',
            "totReads", "nlibraries", "precision"]

    with open(output_file, mode="w", newline='') as csv_out:
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
        cmd='samtools faidx tmp.hairpin.fa'
        run(cmd)

        for x in tqdm.tqdm(mirna_dict,desc="strucVis",disable=None):
            if mirna_dict[x][11]!="No hairpin sequence detected" and mirna_dict[x][11]!="Hairpin structure invalid":
                #Create strucVis plots for each miRNA
                if "Fail" in mirna_dict[x][10]:
                    hp_len=str(mirna_dict[x][9])
                    cmd="strucVis -b alignments/merged.bam -g tmp.hairpin.fa -c "+ x + ":1-" + hp_len + " -s plus -p strucVis/failed/" +x+ ".ps -n " + x
                    run(cmd)
                else:
                    hp_len=str(mirna_dict[x][9])
                    cmd="strucVis -b alignments/merged.bam -g tmp.hairpin.fa -c "+ x + ":1-" + hp_len + " -s plus -p strucVis/passed/" +x+ ".ps -n " + x
                    run(cmd)


    #Create RNAplots depicting miRNA loci
    os.chdir('RNAplots')
    for x in tqdm.tqdm(mirna_dict,desc="RNAplots",disable=None):
        if mirna_dict[x][11] not in fails:
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
    print("Number of submitted candidate miRNAs: " + str(len(mature_dict)))
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


