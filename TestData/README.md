# miRScore Tutorial

## Datasets
There are three data sets provided for testing miRScore.

`TestHairpins.fa/TestmiRNAs.fa dataset` - Dataset containing artificial miRNAs that test each of the different errors that miRScore can encounter when running data.

`candidatemiRNAs.fa/hairpinPrecursors.fa dataset` - Candidate miRNAs acquired through running Mirador on Arabidopsis thaliana sRNA-seq data acquired from PlantSmallRNAgenes.psu.edu.

`miRBase_hairpins.fa/miRBase_miRNAs.fa` - All reported Arabidopsis thaliana miRNAs and corresponding hairpins available on miRBase as of 02/22/2023.


## Instructions

To test these datasets miRScore, install miRScore as directed in the README file in the main repository. 
Download the dataset from the TestData folder that you would like to test miRScore with and place these files in your working directory.

Create a directory called "Fastqs/":
```
mkdir Fastqs
```

Download small sRNA-seq data for Arabidopsis. You can use your own sRNA-seq datasets or download the following:

``` 
cd Fastqs
parallel -j 3 fastq-dump {} ::: SRR5151109 SRR5295847 SRR6192869}
```

Run miRScore in your working directory. For example if you were to run the miRBase dataset, you would use the following command:
```
miRScore --fastq Fastqs/ --mature miRBase_miRNAs.fa --hairpin miRBase_hairpins.fa
```

### Your results will be output into three file sets:

`novel_miRNAS.csv` - miRNAs that miRScore passed  
`failed_miRNAs.csv` - miRNAs that failed at least one criteria for being a miRNA  
`RNAplots` - directory containing hairpin structures with miR/miR* highlight in green (miRNA) and red (miRNA*) 
