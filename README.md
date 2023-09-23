# SeqVerify
SeqVerify is a Python-based command line tool for analysis of whole genome sequencing data for gene-editing verification. It performs insertion site detection, copy number variation (CNV) analysis through CNVPytor, bacterial contamination detection through KRAKEN2 and BRACKEN, and variant calling and filtering aided by SnpEff and SnpSift.

## Install

### Dependencies (& Dependencies of Dependencies)
* BCFtools
* BLAST
* BWA >=0.7
* CNVPytor >=1.3
* HTSLIB
* IDNA
* IGV >= 2.13.2
* Kraken2 >= 2.0
* Matplotlib(-base) >= 2.2
* Numpy
* Python >=3.10
* SAMtools >=1.14
* SciPy
* SNPEff >= 5.1
* SNPSift >= 4.3.1t

### Install through bioconda

```
conda install -c bioconda seqverify
```

## Usage

An example SeqVerify call can be found below:

```
seqverify --output output_name --reads_1 sample_1.fastq --reads_2 sample_2.fastq --genome genome.fa --inexact transgenes.fa --exact commands.txt --database db8gb
```

### Anatomy of a SeqVerify call

SeqVerify has the following standard arguments:
* ```--output``` (Type: String) Used as the identifying name for all the files and folders to do with the particular call. E.g. ```--output Sample1``` will cause the output folder to be named "seqverify_Sample1"
* ```--reads_1``` and ```--reads_2``` (Type: String/Path) The paired-read FASTA/FASTQ, or gzipped FASTA/FASTQ source files for the reads. Also accepts paths to the files if they're not in the working folder.
* ```--genome``` (Type: String/Path) Name or path of FASTA file to be used as reference genome for the reads (e.g. [CHM13](https://github.com/marbl/CHM13#downloads)).
* ```--inexact``` (Type: String/Path) Names or paths of FASTA files containing the sequences of the markers to detect (transgenes, unwanted plasmids, etc.). Accepts more than one if necessary, space-separated. Can also be left blank; not mutually exclusive with ```--exact```. 
* ```--exact``` (Type: String/Path) Name or path to a valid command file for insertion of markers where the insertion site is known. Further details on the construction of a valid command files are given below. Only accepts one command file (but a command file can have multiple commands, so this will not restrict analysis). Can be left blank; not mutually exclusive with ```--inexact```. 

SeqVerify has the following optional arguments:
##### Performance
* ```--threads``` (Type: Integer) Determines how many CPU threads are used in the multithreaded portion of the pipeline. Set to 1 by default.
* ```--max_mem``` (Type: String) Determines what the maximum amount of memory to be used in the indexing of the reference genome should be. Expects an integer followed by "M" or "G" (case-sensitive) for "Megabytes" or "Gigabytes" respectively. Set by default to "16G", 16 Gigabytes.
* ```--start``` (Type: String) Determines where the pipeline starts. Can be set to "beginning" (the default option, runs the entire pipeline), "align" (skips the creation of the augmented genome and the output folders to start at the alignment process), "markers" (skips to the creation of the insertion site readout), "cnv" (skips to the CNV analysis), "plots" (skips to the generation of the CNV plots), "kraken" (skips to the microbial contamination analysis), "variant" (skips to SNV analysis). This option may be useful for core optimization on clusters: e.g. "beginning" is a single-core operation, "align" is a multi-core operation, so on a job scheduler a user could set a job dependent on the other and only use multiple cores when rquired. It may also be used to avoid having to restart the pipeline from scratch should the hardware or the software fail for any reason, or to perform further analysis on data that SeqVerify has already process (e.g. add KRAKEN2 analysis when it wasn't initially requested).
##### KRAKEN2
* ```--kraken``` Determines if KRAKEN2 analysis is to be performed on the sample or not. If set, requires ```--database``` option in order to work correctly.
* ```--database``` (Type: String/Path) Path to valid KRAKEN2 database, or, if KRAKEN2 environmental variables are set, the name of the database. Only needed if ```--kraken``` set. No default setting.
##### Insertion Site Detection
* ```--granularity``` (Type: Integer) Determines how large (in bp) insertion site bins are, i.e. how far apart two insertions can be in order to count as the same insertion site. Set by default to 500, a value of 1 means all insertions on different coordinates will be counted as different insertion sites and show up separately on the readout.
* ```--min_matches``` (Type: Integer) Determines the minimum number of matches/insertions required for an insertion site to appear on the readout. Set by default to 1, i.e. shows any insertion, some of which may be false positives due to repetitive DNA or similar variables.
##### CNVPytor
* ```--bin_size``` (Type: Integer) Determines how many bp are binned together for the purposes of copy number variation detection. Default is 100000 (resulting in 100kbp bins), the minimum value is 1, but anything below the original read length (150 in the test data) will yield meaningless data.
* ```--manual_plots``` Turns off IGV screenshots for coverage plots of transgenes, uses an internal matplotlib script instead. 
##### Variant Calling
* ```--variant_calling``` Turns on the variant calling portion of the pipeline, and takes three space-separated arguments: the name/path of a genome compatible with the VCF annotation DB used (e.g. hg38 for ClinVar), the path to the snpEff config file, and the name/path to a valid annotation DB (e.g. ClinVar)
* ```--variant_intensity``` (Type: String) Includes all variants at or above a certain severity level in the final variant readout. Can be set to (lowest) "MODIFIER", "LOW", "MODERATE", or "HIGH" (highest). Set to "MODERATE" by default.
##### Other
* ```--keep_temp``` If enabled, keeps the temporary files created during the pipeline's execution. It is off by default. We do not recommend turning it on since temp files can reach upwards of 50-100GB depending on the read coverage.
* ```--download_defaults``` If enabled, downloads the default genomes and databases to the working directory: [T2T-CHM13v2.0](https://github.com/marbl/CHM13#analysis-set), intended to be used in ```--genome```, [GRCh38/hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/) for variant calling, and the [8GB PlusPFP](https://benlangmead.github.io/aws-indexes/k2) KRAKEN2 database (placed in a new folder named seqverify_database). Kills the program after downloading these, so anything after ```seqverify --download_defaults``` is ignored.

## Output

SeqVerify will output two folders, ```seqverify_output``` and ```seqverify_temp_output``` where "output" is the variable set in the ```--output``` argument. 

```seqverify_temp_output``` contains files that can be deleted after the pipeline is run, that are generated and used by the pipeline (SAM headers, .BED files, etc.).

```seqverify_output``` contains the following files:

```seqverify_output_markers.bam``` and ```seqverify_output_markers.bam.bai``` A BAM file containing the reads given in ```--reads_1``` and ```--reads_2``` realigned to the reference genome augmented with the marker genomes and its corresponding index.

```seqverify_readout.txt``` The human-readable readout resulting from the insertion site detection portion of the pipeline. The outermost indentation represents the marker, the second indentation is the human chromosome it aligns to, and the third layer denotes the in-chromosome coordinate of the insertion site, and how many times the marker was aligned to it (with the minimum value being specified in ```--min_matches```).

```output.pytor``` The binaries output by CNVPytor. Can be used to generate further plots of specific regions if needed (refer to CNVPytor docs)

```output.global.0000.png``` The Manhattan plot of copy number across the sample as output from ```output.pytor```.

```fig_marker.png``` A copy number histogram of the transgene in the genome. Default bin size is 30bp, and one histogram is output for every sequence specified in ```--markers```. 

If KRAKEN2 analysis is being performed, the pipeline will also output the following files in ```seqverify_output```:

```classified_seqs_output.kreport``` A human-readable report of the microbial sequences detected by KRAKEN2 in text format. Typically species with just a few reads (<10) can be ignored, but the presence of more than this could indicate contamination.

```classified_output.kraken``` The KRAKEN2 binary output files used to generate the report.

```classified_seqs_output_1.fq``` and ```classified_seqs_output_2.fq```, FASTQ files containing the sequences that were classified in the KRAKEN2 database.


If SNV Analysis is being performed, the pipeline will output the following files:

```seqverify_output.ann.vcf``` is the annotated VCF output of all variants present in the sample reads provided. 

 ```seqverify_output_variants.tsv``` is the filtered readout prepared by SeqVerify containing information about all variants above a set severity level, along with any loss of function or homozgosity data.

## Example input and output

Example of WGS data with copy number abnormalities at the beginning of chromosome 12 (CNVPytor Manhattan plot):

![S09_CNV](https://github.com/mpiersonsmela/SeqVerify/assets/20324516/8394f684-c384-4e44-9758-5019443747dd)

Example of the KRAKEN report showing Mycoplasm contamination in the same sample:

![contamination](https://github.com/mpiersonsmela/SeqVerify/assets/20324516/a66666d5-8c4a-4fb8-86b3-3811ffcebb46)

Example of insertion site detection, with tdTomato being found in chromosome 5 of a different sample, and possibly other transgenes (even though given the low number of matches it is more likely they are false positives due to repetitive DNA):

![insertion](https://github.com/mpiersonsmela/SeqVerify/assets/20324516/5010a769-ca21-4b59-a3c2-bfab12b6b234)

## Tutorial

After installing SeqVerify from Bioconda, you will need to assemble the following files to run it:
* Genome sequencing reads, in separated forward and backward FASTQ/FASTA files. SAM format files can be separated using a command like ```samtools fastq``` and then used for the pipeline. These will be referred to as ```r1.fq``` and ```r2.fq```.
* Reference genome, in FASTA format: the pipeline was mainly developed for [CHM13](https://github.com/marbl/CHM13#downloads), but older genomes such as [HG38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) were also tested with good results. This will be referred to as ```genome.fa```.
* Transgene/Marker genome(s), in FASTA format. This tutorial will use two transgenes (even though any number can be used), which will be referred to as ```transgene1.fa``` and ```transgene2.fa```.
* (If KRAKEN is enabled) A valid KRAKEN2 database. More information about how to make a valid KRAKEN2 database can be found [here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases). 

You should then pick a directory to run the pipeline in, from which you'll call the ```seqverify``` command. This will create two folders inside of it, the output folder and the temp folder (the latter of which will be deleted by default). 

```seqverify --output output_name --reads_1 r1.fastq --reads_2 r2.fastq --genome genome.fa --marker_sources transgene1.fa transgene2.fa (--kraken True --database database)```

The pipeline will then run: this will usually take a few hours, with the most time-intensive steps being BWA alignment and genome indexing. The output files as described above will be output in the output folder, and everything else in the temp folder will be deleted unless ```--del_temp``` is set to ```F```.


#### Command text file generation

To use the ```--exact``` flag for known insertion sites of transposons, SeqVerify requires the construction of a commands file. The commands file should be a text file and must be formatted as follows:

```CHR:START-END  SEQUENCE```

where the two columns are tab-separated. Every line is a "command" and specifies a modification to the genome. 

CHR is the name of the chromosome that the modification will happen at (e.g. chr1).

START-END is an interval of coordinates to be deleted, not including START.  

SEQUENCE is the sequence to be inserted from START onwards.

For example, the command ```chr5:56644830-56644850	GGCTCTGGCGAGGGCAGAGGAAGTC``` deletes bases 56644831 to 56644850 in chromosome 5 and inserts the sequence ```GGCTCTGGCGAGGGCAGAGGAAGTC``` in their place. 

Pure insertions (i.e. inserting while deleting nothing) can be achieved by putting the same coordinate for both start and end: ```chr1:1-1 AGCT``` deletes nothing and inserts ```AGCT``` after the first base.
Pure deletions (i.e. deletions with no insertions) can be achieved by leaving the sequence field blank: ```chr2:0-10  ``` deletes the first 10 bases in chromosome 2 and does not replace them with anything.

Multiple commands are allowed in one command file, and seqverify automatically handles interactions between commands on the same chromosome (i.e. a command's coordinates changing because of a previous command), so all the end user needs to do is use the coordinates straight from their source without any adjustment or calculation. 

```seqverify --output output_name --reads_1 r1.fastq --reads_2 r2.fastq --genome genome.fa --inexact transgene1.fa transgene2.fa (--kraken --database database) --threads T --start 2```

Where ```T``` is the number of threads available for the multi-threaded portion.


