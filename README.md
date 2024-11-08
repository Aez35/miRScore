# miRScore

# Description

miRScore is a miRNA validation tool developed to analyze novel miRNAs prior to submission to [miRBase](https://www.mirbase.org/). This tool can be used to determine whether miRNA loci are of high confidence based on the criteria outlined in [Axtell and Meyers (2018)](https://pubmed.ncbi.nlm.nih.gov/29343505/). Users submit mature miRNA, precursor sequences, and sRNA-seq data to miRScore. miRScore will then score each miRNA locus and output several files analyzing the submitted loci. The primary results can be found in the miRScore_Results.csv file, which contains positional and read information, as well as a pass/fail result for each miRNA locus.

### **Determining what is a 'miRNA'**

What constitutes a miRNA is based on several structural and quantitative characteristics. For plants, this criteria can be found in [Axtell and Meyers (2018)](https://pubmed.ncbi.nlm.nih.gov/29343505/). For animals, miRNA criteria is based on those outlined in [Rolle et al. (2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4816427/). Here is a short summary of the criteria outlined in the manuscript:

|microRNA Criteria                                                                |
|:----------------------------------------------------------------------:|
|Less than 5 mismatches in miRNA duplex                                  |
|There must not be more than 3 nucleotides in asymmetric bulges          |
|2 nucleotide 3' overhang of miRNA duplex                                |
|No secondary stems or large loops in miRNA duplex                       |
|Hairpin sequence should be less than 200(animals)/300(plants) nts       |
|miRNA must be between 20-24 nucleotides for plants or 19-25 for animals |
|Confirmation of the miRNA and miRNA* by small RNA-sequencing            |
|At least 75% of reads that map to hairpin must come from miRNA duplex   |
|At least 10 reads within an individual library must map to the miRNA duplex   |

miRScore scores miRNAs based on these conditions, but with some leniency for certain criteria. These exceptions include *miRNAs between 23/24 nt and those with precursors greater than 300 nt*. Instead of treating these as a failed locus, miRScore will report a flag (see flags section for details) indicating that miRNA locus contains one or both of these characteristics. If all other criteria are met, that locus will pass. In addition to these leniencies, miRScore utilizes a read floor instead of requiring multiple libraries have sufficient reads of a miRNA locus. This means that if at least 10 reads map to the miRNA locus in at least one library, that miRNA will have met the criteria. If all other criteria are met, that miRNA will pass.

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
`cutadapt`>=4.8 https://cutadapt.readthedocs.io/en/stable/guide.html
`strucVis`>=0.8
### Create Conda environment on Mac/Linux environment

```
conda create --name miRScore python viennarna pysam bowtie pandas biopython samtools tqdm strucvis cutadapt
```


## miRScore installation

miRScore is run through the python script `miRScore`. To install miRScore, download the file `miRScore` from github page https://github.com/Aez35/miRScore.

miRScore can be made executable by doing the following:

Place the `miRScore` script in your conda environment. To locate this loaction, use the following command:
```
conda info --envs
```

You may make `miRScore` executable and place in the bin directory of your miRScore conda environment:
```
chmod +x miRScore
mv miRScore /Users/zvaricka/miniconda3/envs/miRScore/bin
```

## Required

miRScore requires two FASTA files and multiple small RNA-seq libraries in order to run.   

* `-mature maturefile`: FASTA file containing mature miRNA sequences of proposed novel miRNAs for scoring.   
* `-hairpin hairpinfile`:  FASTA file containing the hairpin sequences in which the mature miRNAs can be found.  
* `-fastq fastqfiles`:  Small RNA-seq libraries in FASTQ format. If libraries are untrimmed, you must use option '-autotrim' to trim libraries. You may supply a key for trimming using option '-trimkey'. Default trimkey is ATH-miR166a for plants and HSA-let-7a for animals. Ideally, all FASTQ files should be **unique individual libraries**. If merged libraries are provided, be sure not to include any individual libraries present in the merged file, as this may cause issues with read counting. It is recommended to avoid the use of merged libraries when possible to ensure highest confidence in miRNA validation. See FAQ for more details.

**Please note: the sequence identifier of each miRNA (i.e. ATH-miR173) must match the sequence identifier of the corresponding hairpin precursors in the hairpin file (i.e. ATH-MIR173)! This is not case sensitive.**

## Usage
```
miRScore [-help] -fastq FASTQFILES.fq/fastq -mature MATUREFILE.fa -hairpin HAIRPINFILE.fa [-star STARFILE.fa] [-mm] [-n NAME] [-autotrim] [-trimkey KEY] [-threads THREADS] -kingdom plant/animal [-out outputDirectory, default= miRScore_output/]
```

## Options

|Option     |Description                                                     |
|:---------:|:-----------------------------------------------------------:   |
|help       | prints help message                                            |
|mature     | FASTA file containing mature and star sequences of novel miRNAs|
|hairpin    | FASTA file containing hairpin precursor sequences of mirnas    |
|fastq      | List of small RNA-seq libraries in FASTQ format                |
|mm         | Allow up to 1 mismatch in miRNA reads                          |
|n          | Specify a name to be added at the beginning of output files    |
|threads    | Specify number of threads to use during trimming step          |
|autotrim   | Trim fastq files using cutadapt                                |
|trimkey    | miRNA sequence for detecting adapter sequences for trimming    |
|kingdom    | Specify either 'plant' or 'animal'                             |
|out        | Specify output directory. Defaults to 'miRScore_output/'       |


### Bowtie
Bowtie has been configured to run the following options.

- a, Report all valid alignments
- v, Allow 0 or 1 mismatch depending on option 'mm'.

Rather than mapping to a genome, miRScore uses Bowtie to map each read to the MIRNA hairpin.


# Example

For this example, miRScore was tested on known _Arabidopsis thaliana_ miRNAs downloaded from miRBase and small RNA-sequencing data downloaded from [PlantSmallRNAgenes](https://plantsmallrnagenes.science.psu.edu/genomes.php?id=50).

1. Download the two FASTA files provided in [TestData](https://github.com/Aez35/miRScore/tree/main/TestData) folder: `ath_miRBase_miRNAs.fa` and `ath_miRBase_precursors.fa`. Place these in your working directory.

2. Create a directory for small RNA-seq libraries and retrieve FASTQ files.

```
mkdir fastqs
cd fastqs
fasterq-dump SRR3222443 SRR3222444 SRR3222445
cd ..
```
3. Run miRScore
```
miRScore -mature ath_miRBase_miRNAs.fa -hairpin ath_miRBase_precursors.fa -fastq fastqs/* -kingdom plant -autotrim
```

# Output
miRScore has six outputs:

* `miRScore_Results.csv`: A csv file containing the results for all proposed miRNAs. Do NOT edit if submitting to miRBase.
* `reads.csv`: A csv file containing read counts for miR, miR*, and total reads mapped to hairpin in each submitted library.
* `alt_miRScore_Results.csv`: A csv file containing the results for alternative miRNA suggestions on a failed miRNA hairpin.
* `alt_reads.csv`: A csv file containing read counts for alternative miR, miR*, and total reads mapped to hairpin in each submitted library
* `RNAplots`: Directory containing PDF images for each miRNA secondary structure including miR(red)/miR*(green) annotation.
* `strucVis`: Directory containing PDF images depicting miRNA secondary structure and small RNA read depth for each miRNA locus.

## miRScore_Results.csv

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
* `totReads` - Total reads mapped to hairpin precursor across all libraries. If a miRNA passes, totReads will report total read counts for only the libraries the miRNA passed in. If the miRNA failed, totReads will report the summed total reads across all submitted libraries. Please note this does not always equal mReads + msReads. This is ALL reads that map to the hairpin.
* `numLibraries` - Number of libraries that reads were report from for mReads, msReads, and totReads.
* `precision` - Measure of precision of how many reads map to the miRNA duplex compared to other positions in the precursor. Found by totaling the number of reads for miR and miR*, then dividing by the sum of all reads mapped to the hairpin precursor.  *(mReads + msReads) / totReads*

In addition to criteria outlined in Axtell and Meyers (2018), miRScore uses a read floor of 10 reads within a single library when evaluation all miRNAs, but no longer requires replication. This means that at least 10 combined miR/miR* reads must be detected in a single library in order for a miRNA to 'Pass'. At least one miR and miR* must still be detected within a library. For example, 9 miR reads and 2 miR* reads in a library will meet this criteria, while 10 miR reads and 0 miR* reads will not.

### **Flags**
Flags are reported to help the user determine what criteria a MIRNA locus did not meet. Most flags reported lead to a hard fail (i.e. 'Precision less than 75%', 'Less than 10 reads'). Two flags ('23/24 nt' and 'Precursor > 300 nt') do not cause a MIRNA locus to fail, and instead are reported to indicate to the user that the miRNA or the miRNA precursor have longer lengths and should be evaluted carefully.

#### List of potential flags
* `hairpin structure invalid` - The miRNA hairpin secondary structure does not allow for the indexing of the miR duplex. This may be due to a large bulge or the miRNA is found too close to the start of the hairpin, not allowing for a 2nt 3' overhang of the miR*.
* `mature structure invalid` - The mature miRNA is likely within a terminal loop structure. This causes the mature miRNA to fold back on itself and fail to meet the criteria.
* `More than 5 mismatches` - More than 5 nucleotides are mismatched between the miRNA and miRNA* sequences.
* `23/24 nt miRNA` - The miRNA/miRNA* is 23 or 24 nucleotides in length.
* `Sequence contained characters besides U,A, T, G, or C` - One of the user-provided sequences contained letters besides A, T, G, C, or U.
* `asymmetric bulge greater than 3` - There is an asymmetric bulge greater than 3 nucleotides in the miR duplex.
* `miRNA not found in hairpin sequence` -The miRNA was not detected in the precursor sequence provided by the user.
* `hairpin is less than 50 basepairs` - The user-provided hairpin sequence is less than 50 basepairs.
* `miRNA multimaps to hairpin` - The miRNA or miRNA* provided by the user multimaps to hairpin sequence.
* `No miR duplex reads detected` - Reads were not detected in the libraries provided for either the miR, miR*, or both.
* `Less than 10 reads` - The miRNA locus had less than 10 combined miR/miR* reads in a single library. The miRNA does not meet the read floor.
* `Precision less than 75%` - The precision (miR* reads + miR* reads / total reads mapped to hairpin) did not meet the required 75% in a single library.
* `No hairpin detected` -The hairpin of the miRNA could not be found in the precursor file provided by the user.
* `precursor > 200/300 nt` - The user provided miRNA precursor has a length greater than 200/300 nucleotides.
* `No 2nt 3' overhang` - The user-provided star sequence did not meet the criteria for a 2nt overhang on the 3' end of the miRNA duplex.
* `NA` - The miRNA meets all miRScore criteria.

## reads.csv
A csv file containing the reads reported in each library of the miRNA, miRNA*, and total reads mapped to each hairpin. Also includes a Pass/Fail result and flag column. A 'Pass' indicates the MIRNA locus that meets the read floor (greater than 10 reads map to MIRNA locus), has at least one read mapped to both miR and miR*, and has a precision greater than 75% in the reported library. A 'Fail' indcates one of these criteria are not met.

## alt_mirna_results.csv
When a miRNA 'fails' miRScore, meaning one or more of the criteria are not met, that locus will be reavaluated based on sequencing data. This revaluation is done using the reads mapped to the hairpin associated with the failed miRNA. miRScore will use samtools to select the most abundant read in the alignments/merged.bam file with a sequence between 20-24 nt that mapped to the failed miRNA hairpin. miRScore will then score that sequence as an alternative miRNA. Any of these alternative miRNAs that pass will be reported in alt_miRScore_Results.csv. These alternative miRNAs are not intended to be final confirmation of a miRNA locus, and should be evaluated with care.
## alt_reads.csv
A csv file containing the reads reported in each library of the alternative miRNA, alternative miRNA*, and total reads mapped to each hairpin. Also includes a Pass/Fail result. A 'Pass' indicates the alternative miRNA meets the read floor and has a precision greater than 75% in that library. A 'Fail' indcates one of these criteria are not met.

## RNAplots

Directory containing .ps files of miRNA secondary structure. This structure is predicted using Vienna RNAfold. For more information on RNA predicted structures, please visit [ViennaRNA page](https://www.tbi.univie.ac.at/RNA/index.html).

Each file contains the hairpin sequence for a miRNA, as well as annotation for the miRNA and predicted miRNA* sequences. The miRNA is outlined in red, while the miRNA* is outlined in green.  
<img width="344" alt="image" src="https://github.com/Aez35/miRScore/assets/88506527/a9b8fa80-084c-4716-b4a4-11189d34bbf5">

## strucvis_plots

Directory containing post-script images generated by strucVis of miRNA secondary structure including small RNA depth of coverage.  
<img width="623" alt="image" src="https://user-images.githubusercontent.com/88506527/230966103-5067eabc-1f5d-4625-b458-9e9fff51d05f.png">


# FAQ
### 1. How does miRScore handle reporting reads?
miRScore reports read counts for the miRNA duplex, totReads, and precision based on which sRNA-seq libraries pass or fail. If a user submits 10 libraries, and all the expression criteria are met in only 2, `mreads`, `msReads`, `totReads`, and `precision` will reflect only those 2 libraries, while the remaining 8 will not be included. Read counts for each locus in each of the 10 libraries can still be found in 'reads.csv' file. If either expression criterion is not met in any individual library, the miRNA will Fail and miRScore will report these values for all 10 libraries.

### 2. Why do some of my miRNAs have 'NA's in certain columns?
miRNA secondary structure plays a big role in determining many of the miRNA criteria. If the hairpin structure forms secondary stem loops, large bulges, or if the miRNA indexes to the stem loop or large loop, the miRNA duplex can be impossible to determine. Without this duplex, it's not possible to score a miRNA locus. In this instance, miRScore will report only the user-provided information and report 'NA's.

### 3. Why is my miRNA failing due to no reads detected, but I see reads reported in both 'mReads' and 'msReads'?
miRscore counts reads in **individual sRNA-seq libraries**. In order for a miRNA locus to pass, there must be reads mapping to both the miRNA and miRNA* positions **within a single library**. This is to prevent false positive validation of miRNAs with a single read in dozens of libraries that may occur if a merged library was used.

miRScore reports three read count values in the results file: `mReads`, `msReads`, and `totReads`. `mReads` and `msReads` are counted using the reported locations of the miRNA/miRNA* on the hairpin +/- 1 position to account for variance. If no reads are detected at one or both of these positions, miRScore will flag that locus as "No mature or star reads detected". In some cases, there are reads mapping to other regions of the hairpin, but not to the miRNA duplex. If this happens, you may see `mReads` or `msReads` reported as 0 while `totReads` is greater than 0. This is because `totReads` is the total number of reads that map to the hairpin, not just the miRNA duplex.

In other cases, you may see "No mature or star reads detected", but 'mReads' and 'msReads' is not reported as 0. This will only occur if you have multiple libraries submitted and miRNA/miRNA* reads are not detected within a single library for ANY library. Meaning, each library contains only miRNA OR miRNA* reads, but not both. If all libraries fail expression criteira, miRScore will report combined read counts for all of the failed libraries. To inspect the read counts in individual libraries, please see the 'reads.csv' file output by miRScore.


### 4. Why is the use of merged libraries not suggested?
As stated in FAQ 3, this is to prevent false positive validation of a MIRNA locus with a single read in dozens of libraries that may occur if a merged library was used. 

The use of merged libraries is not suggested because one of the primary criteria miRScore uses as validation of miRNAs is ensuring both the miRNA and miRNA* is expressed within a single library. This provides a level of confidence that the miRNA is expressed at detectable levels within the sample, and that the reads are in fact miRNA reads with the characteristic two-stack mapping pattern. In some instances, miRNAs annotated by existing sRNA annotation tools only pass when the merged library is used, as it is common practice to use merged libraries to annotate miRNAs. Please be assured that miRScore will still handle merged libraries accurately, but the miRNA and miRNA* abundance reported in the results will reflect the overall abundance and not that of individual libraries. Please see FAQ 1 on how reads are reported when multiple libraries are submitted.




