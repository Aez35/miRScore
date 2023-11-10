# miRScore

# Description

miRScore is a miRNA validation tool developed to analyze novel plant miRNAs prior to submission to [miRBase](https://www.mirbase.org/). This tool can be used to determine whether plant miRNA loci are of high confidence based on the criteria outlined in [Axtell and Meyers (2018)](https://pubmed.ncbi.nlm.nih.gov/29343505/). Users submit mature miRNA and precursor sequences as well as small RNA-seq data to miRScore, which will score each miRNA locus and output a results.csv file containing information and a pass/fail result for each locus.

**Determining what is a 'miRNA'**

Plant miRNAs are detemined based on several structural and quantitative characteristics. This criteria can be found in [Axtell and Meyers (2018)](https://pubmed.ncbi.nlm.nih.gov/29343505/). Here is a short summary of criteria that miRScore uses to determine if a submitted sequence is a miRNA:

|Criteria                                                                |
|:----------------------------------------------------------------------:|
|Less than 5 mismatches in miRNA duplex                                  |
|There must not be more than 3 nucleotides in asymmetric bulges          |
|2 nucleotide 3' overhang of miRNA duplex                                |
|No secondary stems or large loops in miRNA duplex                       |
|miRNA precursor sequence must be less than 300 nucleotides              |
|miRNA must be between 20 to 24 nucleotides                              |
|Confirmation of the miRNA and miRNA* by sRNA-seq only                   |
|At least 75% of reads that map to precursor must come from miRNA duplex |

miRScore scores miRNAs based on these conditions, but with some leniency for certain criteria. These exceptions include *miRNAs between 23/24 nt and those with precursors greater than 300 nt*. Instead of treating these as a failed locus, miRScore will report a flag (see flags section for details) indicating that miRNA locus contains one or both of these characteristics. If all other criteria is met, that locus will pass. In addition to these leniencies, miRScore utilizes a read floor instead of requiring multiple libraries have sufficient reads of a miRNA locus. This means that if at least 10 reads map to the miRNA locus in at least one library, that miRNA will have met the criteria. If all other criteria are met, that miRNA will pass.

# Installation

## Dependencies
All dependencies can be install through `conda`. To install and set up `conda`, follow the instructions at https://bioconda.github.io

`python` >=3.6 https://www.python.org  
`biopython` https://biopython.org/  
`bowtie` >=1.3.1 https://bowtie-bio.sourceforge.net/index.shtml  
`samtools` >=1.16 https://www.htslib.org/  
`viennarna` >=2.5.1 https://www.tbi.univie.ac.at/RNA/documentation.html  
`pysam` >=0.21.0 https://pysam.readthedocs.io/en/latest/api.html  
`tqdm` >=4.65 https://pypi.org/project/tqdm/  
`pandas`>=2.0.0 https://pandas.pydata.org/  

### Create Conda environment on Mac/Linux environment

```
conda create --name miRScore python viennarna pysam bowtie pandas biopython samtools tqdm
```


## miRScore installation

miRScore is run through the python script `miRScore`. To install miRScore, download the file `miRScore` from github page https://github.com/Aez35/miRScore.

miRScore can be made executable by doing the following:

Place the `miRScore` script in your conda environment. To locate this loaction, use the following command:
```
conda info --envs
```

Make `miRScore` executable and place in the bin directory of your miRScore conda environment:
```
chmod +x miRScore
mv miRScore /Users/zvaricka/miniconda3/envs/miRScore/bin
```

## Required

miRScore requires two FASTA files and multiple small RNA-seq libraries in order to run.   

* `-mirnas mirnafile`: FASTA file containing mature miRNA sequences of proposed novel miRNAs for scoring.   
* `-precursors precursorfile`:  FASTA file containing the hairpin sequences in which the mature miRNAs can be found.  
* `(-fastq fastqfiles | -bamfile bamfiles)`: Small RNA-seq libraries, either in FASTQ or BAM format. The specified files should be **unique individual libraries**. If merged libraries are provided, be sure not to include any individual libraries present in the merged file, as this will cause issues with read counting. It is best to avoid the use of merged libraries when possible.

**Please note that the sequence identifier of each miRNA must match the sequence identifier of the corresponding hairpin precursors.**


## Usage
```
miRScore [-help] ([-fastq FASTQFILES]|[-bamfile BAMFILES]) -mirnas MIRNAFILE -precursors PRECURSORFILE [-mm] [-n NAME] [-threads THREADS]
```

## Options

|Option     |Description                                                  |
|:---------:|:-----------------------------------------------------------:|
|help       | prints help message                                         |
|mirnas     | FASTA file containing mature sequences of novel miRNAs      |
|precursors | FASTA file containing hairpin precursor sequences of mirnas |
|fastq      | List of small RNA-seq libraries in FASTQ format             |
|bam        | List of small RNA-seq libraries in bam format               |
|mm         | Allow up to 1 mismatch in miRNA reads                       |
|n          | Specify a name to be added at the beginning of output files |
|threads    | Specify number of threads to use during bam2fastq step      |

### Bowtie
Bowtie has been configured to run the following options.

- a, Report all valid alignments
- v, Allow 0 or 1 mismatch depending on option 'mm'.

Rather than mapping to a genome, miRScore uses Bowtie to map each read to the miRNA hairpin.


# Example

For this example, miRScore was tested on known _Arabidopsis thaliana_ miRNAs downloaded from miRBase and small RNA-sequencing data downloaded from [PlantSmallRNAgenes](https://plantsmallrnagenes.science.psu.edu/genomes.php?id=50).

1. Download the two FASTA files provided in [TestData](https://github.com/Aez35/miRScore/tree/main/TestData) folder: `ath_miRBase_miRNAs.fa` and `ath_miRBase_precursors.fa`. Place these in your working directory.

2. Create a directory for small RNA-seq libraries and retrieve bam files.

```
mkdir bamfiles
cd bamfiles
wget http://plantsmallrnagenes.science.psu.edu/ath-b10/alignments/SRR3222445_trimmedSingle.bam
wget http://plantsmallrnagenes.science.psu.edu/ath-b10/alignments/SRR3222444_trimmedSingle.bam
wget http://plantsmallrnagenes.science.psu.edu/ath-b10/alignments/SRR3222443_trimmedSingle.bam
cd ..
```
3. Run miRScore
```
miRScore -mirnas ath_miRBase_miRNAs.fa -precursors ath_miRBase_precursors.fa -bamfile bamfiles/*
```

# Output
miRScore has three outputs:

* `miRScore_Results.csv`: A csv file containing the results for all proposed miRNAs. Do NOT edit if submitting to miRBase.
* `reads.csv`: A csv file containing read counts for miR, miR*, and total reads mapped to hairpin in each submitted library.
* `alt_miRScore_Results.csv`: A csv file containing the results for alternative miRNA suggestions on a failed miRNA hairpin.
* `alt_reads.csv`: A csv file containing read counts for alternative miR, miR*, and total reads mapped to hairpin in each submitted library
* `RNAplots`: Directory containing post-script images for each miRNA secondary structure including miR(red)/miR*(green) annotation.
* `strucVis`: Directory containing post-script images depicting miRNA secondary structure and small RNA read depth for each miRNA locus.

### miRScore_Results.csv

This file contains the results for all candidate miRNAs submitted to miRScore. Do NOT edit this file if submitting to miRBase. miRBase submission requires specific formatting, specified by miRScore output. Editing this file may cause issues in submission.

* `name` - miRNA name provided by user.  
* `mSeq` - miRNA sequence.  
* `mLen` - miRNA length.  
* `mStart` - Start position of miRNA on hairpin precursor.  
* `mStop` - Stop position of miRNA on hairpin precursor.   
* `msSeq`- miRStar sequence.  
* `msLen` - miRStar length.  
* `msStart` - Start position of miRStar on hairpin precursor.  
* `msStop` - Stop position of miRStar on hairpin precursor.  
* `precSeq` - Hairpin precursor sequence.  
* `precLen` - Length of hairpin precursor.  
* `result` - Pass/Fail result of miRScore for miRNA.  
* `flags` - Criteria that was not met by the miRNA duplex according to Axtell and Meyers (2018).  
* `mismatch` - Number of mismatches in miRNA duplex.  
* `mReads` - Sum of reads in libraries that mapped to miRNA locus +/- 1 position. If a miRNA passes, mReads will report miR read counts for only the libraries the miRNA passed in. If the miRNA failed, mReads will report the summed miR reads across all submitted libraries.
* `msReads` - Sum of reads in all libraries that mapped to miRStar locus +/- 1 position. If a miRNA passes, msReads will report miR* read counts for only the libraries the miRNA passed in. If the miRNA failed, msReads will report the summed miR* reads across all submitted libraries.
* `totReads` - Total reads mapped to hairpin precursor across all libraries. If a miRNA passes, totReads will report total read counts for only the libraries the miRNA passed in. If the miRNA failed, totReads will report the summed total reads across all submitted libraries.
* `numLibraries` - Number of libraries that reads were report from for mReads, msReads, and totReads.
* `precisionScore` - Measure of read precision at miRNA locus. Found by totaling the number of reads for miR and miR*, then dividing by the sum of all reads mapped to the hairpin precursor.  *(mReads + msReads) / totReads*

In addition to criteria outlined in Axtell and Meyers (2018), miRScore uses a read floor of 10 reads within a single library when evaluation all miRNAs, but no longer requires replication. This means that at least 10 combined miR/miR* reads must be detected in a single library in order for a miRNA to 'Pass'. At least one miR and miR* must still be detected within a library. For example, 9 miR reads and 2 miR* reads in a library will meet this criteria, while 10 miR reads and 0 miR* reads will not.

  #### flags


### alt_miRScore_Results.csv
When a miRNA 'fails' miRScore, meaning that not all 

### RNAplots

Directory containing .ps files of miRNA secondary structure. This structure is predicted using Vienna RNAfold. For more information on RNA predicted structures, please visit [ViennaRNA page](https://www.tbi.univie.ac.at/RNA/index.html).

Each file contains the hairpin sequence for a miRNA, as well as annotation for the miRNA and predicted miRNA* sequences. The miRNA is outlined in red, while the miRNA* is outlined in green.  
<img width="623" alt="image" src="https://user-images.githubusercontent.com/88506527/230966219-33cc6f62-a1dc-46fe-8541-176ef43cca88.png">

### strucvis_plots

Directory containing post-script images generated by strucVis of miRNA secondary structure including small RNA depth of coverage.  
<img width="623" alt="image" src="https://user-images.githubusercontent.com/88506527/230966103-5067eabc-1f5d-4625-b458-9e9fff51d05f.png">
