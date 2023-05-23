# miRScore

# Description

miRScore is a miRNA scoring tool that was developed to analyze novel miRNAs prior to submission to [miRBase](https://www.mirbase.org/). This tool can be used to determine whether a miRNA loci is of high confidence based on the criteria outlined in [Axtell and Meyers (2018)](https://pubmed.ncbi.nlm.nih.gov/29343505/). Users submit mature miRNA and precursor sequences as well as small RNA-seq data to miRScore, which will output a .csv file containing information on each miRNA loci and a pass/fail result.

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
* `(-fastq fastqDirectory | -bam bamfilesDirectory)`: Directory containing small RNA-seq libraries, either in FASTQ or BAM format. The specified directory should only contain **unique individual libraries**. miRScore will include any file with a '.bam' or '.fastq' extension according to which option is used. Please be sure to remove any merged library files before running miRScore.

**Please note that the sequence identifier of each miRNA must match the sequence identifier of the corresponding hairpin precursors.**


## Usage
```
miRScore [-help] ([-fastq fastqDirectory]|[-bam bamfilesDirectory]) -mirnas MIRNAFILE -precursors PRECURSORFILE [-mm]
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

For this example, miRScore was tested on known _Arabidopsis thaliana_ miRNAs downloaded from miRBase and small RNA-sequencing data downloaded from [PlantSmallRNAgenes](https://plantsmallrnagenes.science.psu.edu/genomes.php?id=50). Before you begin, download the two test data FASTA files provided in [TestData](https://github.com/Aez35/miRScore/tree/main/TestData) folder: `ath_miRBase_miRNAs.fa` and `ath_miRBase_precursors.fa`. Place these in your working directory.


1. Create a directory for small RNA-seq libraries and retrieve bam files.

```
mkdir bamfiles
cd bamfiles
wget http://plantsmallrnagenes.science.psu.edu/ath-b10/alignments/SRR5295846_trimmedSingle.bam
wget http://plantsmallrnagenes.science.psu.edu/ath-b10/alignments/SRR5151110_trimmedSingle.bam
cd ..
```
2. Run miRScore
```
miRScore -mirnas ath_miRBase_miRNAs.fa -precursors ath_miRBase_precursors.fa -bam bamfiles/
```

# Output
miRScore has three outputs:

* `miRScore_Results.csv`: A csv file containing the results for all proposed miRNAs
* `RNAplots`: Directory containing post-script images for each miRNA secondary structure including miR(red)/miR*(green) annotation.
* `strucVis`: Directory containing post-script images depicting miRNA secondary structure and small RNA read depth for each miRNA loci.

### miRScore_Results.csv

This file contains the results for all candidate miRNAs submitted to miRScore.

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
* `mismatch` - Number of mismatches in miRNA.  
* `mReads` - Sum of reads in all libraries that mapped to miRNA loci +/- 1 position.  
* `msReads` - Sum of reads in all libraries that mapped to miRStar loci +/- 1 position.  
* `totReads` - Total reads mapped to hairpin precursor across all libraries.  
* `precision_score` - Measure of read precision at miRNA loci. Found by totaling the number of reads for miR and miRStar, then dividing by the sum of all reads mapped to the hairpin precursor.  *(mReads + msReads) / totReads*


### RNAplots

Directory containing .ps files of miRNA secondary structure. This structure is predicted using Vienna RNAfold. For more information on RNA predicted structures, please visit [ViennaRNA page](https://www.tbi.univie.ac.at/RNA/index.html).

Each file contains the hairpin sequence for a miRNA, as well as annotation for the miRNA and predicted miRNA* sequences. The miRNA is outlined in red, while the miRNA* is outlined in green.  
<img width="623" alt="image" src="https://user-images.githubusercontent.com/88506527/230966219-33cc6f62-a1dc-46fe-8541-176ef43cca88.png">

### strucvis_plots

Directory containing post-script images generated by strucVis of miRNA secondary structure including small RNA depth of coverage.  
<img width="623" alt="image" src="https://user-images.githubusercontent.com/88506527/230966103-5067eabc-1f5d-4625-b458-9e9fff51d05f.png">
