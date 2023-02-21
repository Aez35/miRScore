# miRScore

# Description

miRScore is a miRNA scoring tool developed to analyze novel miRNAs for submission to miRBase, a miRNA database. The criteria used to determine whether a novel miRNA is of high confidence is outlined in Axtell and Meyers (2018).

miRScore requires two fasta files as input: one of novel miRNAs annotated prior to running miRScore and one of their hairpins. Please note that the identifier of each miRNA must match the identifier of the corresponding hairpin. In addition, this program requires multiple (two or more) fastq files containing evidence of these candidates to determine if a candidate miRNA meets the outlined criteria. 

Running miRScore outputs two csv files:   
`miRNA_candidats.csv` - csv file containing candidates miRNAs that met all criteria and relevant information for each loci (i.e. mature sequence, length of miRNA, etc)   
`failed_miRNAs.csv` - csv file contanins those that failed to meet all criteria and a short description of which criteria it failed on

# Installation

## Dependencies
All dependencies can be install through `conda`. To install and set up `conda`, follow the instructions at https://bioconda.github.io

### Linux or Mac
```
conda create --name miRScore python=3.6 ViennaRNA pysam argparse bowtie pandas biopython regex pathlib
```


## miRScore installation

miRScore is run through the python script `miRScore`. To install miRScore, download the file `miRScore` from github page https://github.com/Aez35/miRScore.

Place the `miRScore` script in your conda environment. To locate this loaction, use the command `conda info --envs`.

For example:
```
chmod +x miRScore
mv miRScore /opt/anaconda3/envs/miRScore
```
    
# Configuration

##Usage
```
miRScore [-h] ([--fastq fastqDirectory]|[--bam bamfileDirectory]) --mature MIRNAFILE --hairpin HAIRPINFILE [-mm]
```

## Options

|Option     |Description                                               |
|:---------:|:--------------------------------------------------------:|
|mature     | Fasta file containing mature sequences of novel miRNAs   |
|hairpin    | Fasta file containing hairpin sequences of novel miRNAs  |
|fastq      | Directory containing two or more fastq files             |
|bam        | Directory containing two or more bam files               |
|mm         | Allow up to 1 mismatch in miRNA reads.                   |

## Bowtie

Bowtie has been configured to run the following options.

- a, Report all valid alignments

Rather than mapping to a genome, miRScore uses Bowtie to map each read to the miRNA precursors.


# Example

Create a directory containing fastq files and run miRSCore.
```
mkdir fastqs
fasterq-dump SRR3222445 -O fastqs/
fasterq-dump SRR3269282 -O fastqs/
python3.6 miRScore.py --mature candidatemiRNAs.fa --hairpin precursors.fa --fastqd fastqs/
```

# Output
miRScore outputs two csv files that contain information about successful and fail candidates.

**novel_miRNAS.csv** contains all miRNAs that had reads present in two or more libraries and who's miR/miR* duplex met requirements of being a miRNA.

**failed_miRNAs.csv** contains the identifier of each miRNA that failed and the reason it failed.
