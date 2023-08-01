# SeqVerify
SeqVerify is a Python-based command line tool for analysis of whole genome sequencing data for gene-editing verification. It performs insertion site detection, copy number variation (CNV) analysis through CNVPytor, and bacterial contamination detection through KRAKEN2 and BRACKEN.

## Install

### Dependencies
* Python >=3.9
* SAMtools >=1.14
* BWA >=0.7
* CNVPytor >=1.3
* Kraken2 >= 2.0
* Matplotlib >= 2.2

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
* ```--start``` (Type: Integer) For core-usage optimization. If set to 1, runs the entire pipeline from the start. If set to 0, only runs the single-threaded parts of the pipeline (i.e. up to and including the indexing of the reference genome, but stopping after that). If set to 2, only runs the multi-threaded parts of the pipeline (i.e. assumes the reference genome is already indexed and continues from there). Set to 1 by default, the logic behind this argument being that, where core optimization may be desirable, the pipeline can be split up into two commands, one that just requests a single thread, and one that requests multiple. Additionally, if one is using the same reference chromosome and markers for multiple samples, it is useless to re-index the genome every time, so setting 0 can be run once for the entire cohort, and setting 2 can be run after that for each individual sample.
##### KRAKEN2
* ```--kraken``` (Type: String/Boolean) Determines if KRAKEN2 analysis is to be performed on the sample or not. Set to False by default, if set to True requires ```--database``` option in order to work correctly.
* ```--database``` (Type: String/Path) Path to valid KRAKEN2 database, or, if KRAKEN2 environmental variables are set, the name of the database. Only needed if ```--kraken``` set to True. No default setting.
##### Insertion Site Detection
* ```--granularity``` (Type: Integer) Determines how large (in bp) insertion site bins are, i.e. how far apart two insertions can be in order to count as the same insertion site. Set by default to 500, a value of 1 means all insertions on different coordinates will be counted as different insertion sites and show up separately on the readout.
* ```--min_matches``` (Type: Integer) Determines the minimum number of matches/insertions required for an insertion site to appear on the readout. Set by default to 1, i.e. shows any insertion, some of which may be false positives due to repetitive DNA or similar variables.
##### CNVPytor
* ```--bin_size``` (Type: Integer) Determines how many bp are binned together for the purposes of copy number variation detection. Default is 100000 (resulting in 100kbp bins), the minimum value is 1, but anything below the original read length (150 in the test data) will yield meaningless data.
##### Other
* ```--del_temp``` If enabled, deletes the temporary files created during the pipeline's execution. It is on by default, and highly recommended, since temp files can reach upwards of 50-100GB depending on the read coverage.
* ```--download_defaults``` If enabled, downloads the default genomes and databases to the working directory: [T2T-CHM13v2.0](https://github.com/marbl/CHM13#analysis-set), intended to be used in ```--genome```, [GRCh38/hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/) for variant calling, and the [8GB PlusPFP](https://benlangmead.github.io/aws-indexes/k2) KRAKEN2 database (placed in a new folder named seqverify_database). Kills the program after downloading these.


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

#### Core optimization

If you are running SeqVerify on a cluster or other powerful machine, you may be interested in optimizing the usage of the cores you use at one. This is where the ```--start``` flag becomes useful: instead of scheduling one ```seqverify``` command as shown above, you can call two, one with ```--start 0``` and a subsequent one with ```--start 2``` to split up the single-threaded and multi-threaded portions of the pipeline. The rest of the arguments should be able to remain the same, such that the two commands you call are:

```seqverify --output output_name --reads_1 r1.fastq --reads_2 r2.fastq --genome genome.fa --marker_sources transgene1.fa transgene2.fa (--kraken True --database database) --start 0```

```seqverify --output output_name --reads_1 r1.fastq --reads_2 r2.fastq --genome genome.fa --marker_sources transgene1.fa transgene2.fa (--kraken True --database database) --threads T --start 2```

Where ```T``` is the number of threads available for the multi-threaded portion.

#### Cohort 

If you are running SeqVerify on multiple samples with the same underlying genome and transgene markers, you can use the ```--start 0``` option combined with ```--del_temp F``` to run the pipeline up to indexing the new genome, to then halt it, manually go into the temp folder, and copy-paste the reference genome and its index into new manually created temp folders for the other samples (making sure the naming convention matches the one described in the main ```seqverify``` file). You will then be able to run ```--start 2``` for all the other samples and skip genome indexing, saving repetitive operations, time, and compute.


