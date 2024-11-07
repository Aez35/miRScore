#!/usr/bin/env python
#miRScore V0.2.1

# This python script was written to address issues in datasets being submitted to miRScore where multiple miRNAs come from a single locus. 
# This occurs with miRBase datasets frequently, where two miRNAs (i.e. osa-miR159a.1 and osa-miR159a.2) have a single hairpin precursor sequence listed (i.e. osa-MIR159a).
# The simple solution is to create a new hairpin precursor FASTA file with a hairpin sequence for each miRNA.
# Running this script will create a "miRScore_adjusted_hairpins.fa" file which will have hairpin precursors that match the offending miRNA sequences.

import argparse, subprocess,sys
from Bio import SeqIO
import os, tqdm, glob
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter

print('')
#Assign arguments for command line
ap = argparse.ArgumentParser()

ap.add_argument("-mature", required = True, help="fasta file of mature miRNA sequence")
ap.add_argument("-hairpin", required = True, help="fasta file of hairpin precursor sequence")
args = ap.parse_args()

#Check that mature and hairpin files exist.
if not os.path.isfile(args.mature):
    msg=(args.mature + " does not exist")
    sys.exit(msg)
if not os.path.exists(args.hairpin):
    msg=(args.hairpin + " does not exist")
    sys.exit(msg)


def check_multiple(hp_dict,mir_dict):
    # Temporary dictionary to hold new entries
    new_entries = {}

    # Loop through each key in hp_dict
    for key in list(hp_dict.keys()):  # Use list(hp_dict.keys()) to avoid modifying during iteration
        # Track if new entries are added for this key
        has_new_entries = False

        # Check for keys in mir_dict that start with the current key from hp_dict
        for sub_key, seq_record in mir_dict.items():

            if sub_key.upper().startswith(key.upper()) and not (sub_key.endswith("-5p") or sub_key.endswith("-3p")):
                # Create a duplicate record for each MIRNA locus with multiple miRNAs.
                new_record = SeqRecord(
                    seq=hp_dict[key].seq,
                    id=seq_record.id.replace("MIR", "miR"),
                    name=seq_record.name.replace("MIR", "miR"),
                    description=seq_record.description.replace("MIR", "miR"),
                    dbxrefs=seq_record.dbxrefs
                )
                # Add the new record to new_entries
                new_entries[sub_key] = new_record
                has_new_entries = True  # Mark that we've added a new entry

        # If new entries were added, remove the original key from hp_dict
        if has_new_entries:
            del hp_dict[key]

    # Update hp_dict with the collected new entries
    hp_dict.update(new_entries)
    if len(new_entries)>0:
        # Specify the output FASTA file path
        output_fasta = "miRScore_adjusted_hairpins.fa"

        # Open the output file and write sequences in single-line FASTA format
        with open(output_fasta, "w") as output_handle:
            for i, record in enumerate(hp_dict.values()):
                # Write the header
                output_handle.write(f">{record.id}\n")
                # Write the sequence without wrapping
                output_handle.write(f"{str(record.seq)}")
                # Add a newline after each record, except the last one
                if i < len(hp_dict) - 1:
                    output_handle.write("\n")

def process_mirnas(input_file, output_file):
    #Convert all U's to T's for mapping.
    subprocess.run(
        f"cat {input_file} | sed -E '/(>)/!s/U/T/g' > {output_file}",
        shell=True
    )
    
def run(cmd) :
#Run subprocess
    proc = subprocess.call(cmd, shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

#__________________________ Code begins ________________________________
def main():
    #Create temporary files converting any U's to T's
    #hairpins
    process_mirnas(args.hairpin, "tmp.hairpin.fa")
    #miRNA
    process_mirnas(args.mature, "tmp.mature.fa")

    #Fetch miRNA precursor sequence matching mature candidate miRNA sequence
    mir_dict1 =  SeqIO.index("tmp.mature.fa", "fasta")
    hp_dict1 = SeqIO.index("tmp.hairpin.fa", "fasta")
    mir_dict = {key: mir_dict1[key] for key in mir_dict1}
    hp_dict = {key: hp_dict1[key] for key in hp_dict1}

    hpLen=len(hp_dict)
    check_multiple(hp_dict,mir_dict)
    newLen=len(hp_dict)

    print("Created new hairpin FASTA file called 'miRScore_adjusted_hairpins.fa' with " + str(newLen-hpLen) + " additional hairpin sequence(s).")
    print('')

    cmd="rm -rf tmp*"
    run(cmd)
    
if __name__ == "__main__":
    main()
