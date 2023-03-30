# miRScore

# Description

miRScore is a miRNA scoring tool that analyzes novel miRNAs. 

*The criteria used to determine whether a novel miRNA is of high confidence is outlined in [Axtell and Meyers (2018)](https://pubmed.ncbi.nlm.nih.gov/29343505/).*

# Installation

## Dependencies
All dependencies can be install through `conda`. To install and set up `conda`, follow the instructions at https://bioconda.github.io

`Python` >=3.6 https://www.python.org  
`Biopython` https://biopython.org/
`Bowtie` >=1.3.1 https://bowtie-bio.sourceforge.net/index.shtml  
`Samtools` >=1.16 https://www.htslib.org/  
`ViennaRNA2` >=2.5.1 https://www.tbi.univie.ac.at/RNA/documentation.html  
`Sratoolkit` https://hpc.nih.gov/apps/sratoolkit.html  
`Pysam` >=0.21.0 https://pysam.readthedocs.io/en/latest/api.html

### Conda installation on Mac/Linux environment

```
conda create --name miRScore python=3.6 ViennaRNA pysam argparse bowtie pandas biopython regex pathlib sra-tools
```


## miRScore installation

miRScore is run through the python script `miRScore.py`. To install miRScore, download the file `miRScore.py` from github page https://github.com/Aez35/miRScore and place the script in your working directory.

    
## Required

miRScore requires two FASTA files and multiple small RNA-seq libraries in order to run.   

* `--mirnas mirnafile`: mature miRNA sequences of proposed novel miRNAs for scoring.   
* `--precursors precursorfile`:  hairpin precursors and one containing the hairpin sequences in which the mature miRNAs can be found.  
* `(--fastq fastqDirectory | --bam bamfilesDirectory)`: Directory containing small RNA-seq libraries, either in FASTQ or BAM format.  

**Please note that the sequence identifier of each miRNA must match the sequence identifier of the corresponding hairpin.**


## Usage
```
python3.6 miRScore.py [--help] ([--fastq fastqDirectory]|[--bam bamfilesDirectory]) --mirnas MIRNAFILE --precursors PRECURSORFILE [-mm]
```

## Options

|Option     |Description                                                  |
|:---------:|:-----------------------------------------------------------:|
|help       | prints help message                                         |
|mirnas     | FASTA file containing mature sequences of novel miRNAs      |
|precursors | FASTA file containing hairpin precursor sequences of mirnas |
|fastq      | Directory containing two or more fastq files                |
|bam        | Directory containing two or more bam files                  |
|mm         | Allow up to 1 mismatch in miRNA reads                       |

### Bowtie

Bowtie has been configured to run the following options.

- a, Report all valid alignments
- v, Allow 0 or 1 mismatch depending on option 'mm'.

Rather than mapping to a genome, miRScore uses Bowtie to map each read to the miRNA hairpin.


# Example

For this example, miRScore was tested on known _Arabidopsis thaliana_ miRNAs downloaded from miRBase. Before you begin, download the `mature.fa` and `hairpin.fa` files from [miRBase](https://www.mirbase.org/ftp.shtml).

1. Extract the *Arabidopsis thaliana* entries using the following commands:
```
grep -A 1 '>ath' mature.fa | grep -v '\-\-' | sed 's/ .*//g'|tr '[:lower:]' '[:upper:]'| sed 's/..P//g' > ath_miRNAs.fa
grep -A 1 '>ath' hairpin.fa | grep -v '\-\-' | sed 's/ .*//g'| tr '[:lower:]' '[:upper:]' > ath_precursors2.fa
```
2. Create a directory for small RNA-seq libraries and retrieve FASTQ files.

```
mkdir fastqs
cd fastqs
fasterq-dump SRR3222443 SRR3222444
cd ..
```
3. Run miRScore
```
python3.6 miRScore.py --mirnas ath_miRNAs.fa --precursors ath_precursors.fa --fastqd fastqs/
```

# Output
miRScore outputs a csv file containing the results for all proposed miRNAs , `miRScore_Results.csv`.


* `Name` - miRNA name provided by user.  
* `mSeq` - miRNA sequence.  
* `mLen` - miRNA length.  
* `mStart` - Start position of miRNA on hairpin precursor.  
* `mStop` - Stop position of miRNA on hairpin precursor.   
* `mStarSeq`- miRStar sequence.  
* `mStarLen` - miRStar length.  
* `msStart` - Start position of miRStar on hairpin precursor.  
* `msStop` - Stop position of miRStar on hairpin precursor.  
* `precSeq` - Hairpin precursor sequence.  
* `precLen` - Length of hairpin precursor.  
* `result` - Pass/Fail result of miRScore for miRNA.  
* `flags` - Criteria that was not met by the miRNA according to Axtell and Meyers (2018).  
* `mismatch` - Number of mismatches on miRNA.  
* `mReads` - Sum of reads in all libraries that mapped to miRNA loci +/- 1 position.  
* `msReads` - Sum of reads in all libraries that mapped to miRStar loci +/- 1 position.  
* `totReads` - Total reads mapped to hairpin precursor across all libraries.  
* `precision` - Measure of read precision at miRNA loci. Found by totaling the number of reads for miR and miRStar, then dividing by the sum of all reads mapped to the hairpin precursor.  *(mReads + msReads) / totReads*


