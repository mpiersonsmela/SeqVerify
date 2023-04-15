#!/usr/bin/env python
#Magnify - A sam parser for insertion site detection
 
import argparse
import os
from seqver_functions import *
from seqver_transgenes import *

parser = argparse.ArgumentParser()
#python magnify_commandline.py --output S04 --reads_1 S04_1.fastq --reads_2 S04_2.fastq --genome_source chm13v2.0.fa --marker_sources T2A-mGreenLantern.fa T2A-tdTomato.fa transposons_in_S04.fa unwanted_plasmids.fa --granularity 500 --threads 15
parser.add_argument('--output', type=str, required=True)
parser.add_argument('--reads_1', type=str, required=True)
parser.add_argument('--reads_2', type=str, required=True)
parser.add_argument('--genome_source', type=str, required=True)
parser.add_argument('--marker_sources', type=str, required=True, nargs='+')
parser.add_argument('--database',type=str,required=True)

parser.add_argument('--granularity',type=int, required=False, default=500)
parser.add_argument('--threads',type=int,required=False, default=1)
parser.add_argument('--max_mem', type=str, required=False, default='16G')
parser.add_argument('--min_matches', type=int, required=False,default=1)
parser.add_argument('--kraken', type=str, required=False,default="F")

parser.add_argument('--bin_size',type=int, required=False, default=100000)

args = parser.parse_args()

output_name = args.output

#The pathFinder function detects whether arguments are paths or names and acts accordingly 
#such that it can handle both cases for ease of use.
reads_1_path = pathFinder(args.reads_1)
reads_2_path = pathFinder(args.reads_2)
genome_source_path = pathFinder(args.genome_source)
marker_sources_path = [pathFinder(i) for i in args.marker_sources]

granularity_threshold = args.granularity
database = pathFinder(args.database)

min_matches = args.min_matches
bin_size = args.bin_size
threads = args.threads

if args.bin_size != 100000:
    if args.bin_size % 100 == 0:
        bin_size = args.bin_size
    else:
        raise ValueError("Custom bin size not divisible by 100")

#Deals with memory logic, converts max_mem argument into number of bytes of memory requested in order to build optimal BWA indexing block size
if args.max_mem[-1] == "M":
    max_mem = int(args.max_mem[:-1])*1000000
elif args.max_mem[-1] == "G":
    max_mem = int(args.max_mem[:-1])*1000000000
elif type(args.max_mem) == 'int':
    max_mem = int(args.max_mem)

base_block = max_mem / 8

#Preliminary Work - Generates the output folder given a specified --output string
#                   as well as the temp folder for the pipeline to work within
os.system(f'mkdir magnify_{output_name}')
os.system(f'mkdir magnify_temp_{output_name}')

#Prep Collated Genome:

#Collates all markers/transgenes into a single FASTA file for ease of handling
for file in marker_sources_path: 
    os.system(f'echo "$(cat {file})" >> magnify_temp_{output_name}/magnify_{output_name}_transgene_list.fa') 

#Copies the reference genome into a new FASTA file that the transgenes will also be added to
os.system(f'cp {genome_source_path} magnify_temp_{output_name}/magnify_{output_name}_collated.fa')

#Adds the transgenes to the new "collated" FASTA file
os.system(f'echo "$(cat magnify_temp_{output_name}/magnify_{output_name}_transgene_list.fa)" >> magnify_temp_{output_name}/magnify_{output_name}_collated.fa')

#Builds the indices for the genome that has now been augmented with the transgenes, necessary for alignment
os.system(f'bwa index -b {base_block} magnify_temp_{output_name}/magnify_{output_name}_collated.fa')


#Aligns the new genome using BWA-MEM
os.system(f'bwa mem -M -t {threads} magnify_temp_{output_name}/magnify_{output_name}_collated.fa {reads_1_path} {reads_2_path} > magnify_temp_{output_name}/magnify_{output_name}_markers.sam')

#Saves all names of transgenes/markers we want to find insertion sites of into a list for later usage and makes it into a string
with open(f'magnify_temp_{output_name}/magnify_{output_name}_transgene_list.fa') as f:
    marker_list = [i[1:].strip() for i in f if i.startswith('>')]
markers = "|".join(marker_list)

#Uses the "markers" string as a regular expression to find reads in the newly-aligned SAM file that align to a transgene and have a mate on a human chromosome
#Exports these reads into a new SAM file, "diff_chr"
os.system(f"samtools view -h -@ {threads} magnify_temp_{output_name}/magnify_{output_name}_markers.sam |gawk '{{if ((($3 ~ /^({markers})/)&&($7 ~ /chr([0-9]+|[XYM])/))||($1 ~ /^@/)) print $0}};' > magnify_temp_{output_name}/magnify_{output_name}_markers_diff_chr.sam")

#Functions imported from magnify_functions, refer to comments there for details on their functioning.
#Takes "diff_chr", generates a .txt human-readable readout of the insertion sites, places it into the output folder.
data = group(f'magnify_temp_{output_name}/magnify_{output_name}_markers_diff_chr.sam')
insertions = compress(data,granularity_threshold)
readout(insertions,marker_list,min_matches)
os.system(f"mv magnify_readout.txt magnify_{output_name}")

os.system(f"samtools sort -@ {threads} magnify_temp_{output_name}/magnify_{output_name}_markers.sam > magnify_{output_name}/magnify_{output_name}_markers.bam")
os.system(f"samtools index -@ {threads} magnify_{output_name}/magnify_{output_name}_markers.bam")

os.system(f'cnvpytor -root magnify_{output_name}/magnify.pytor -rd magnify_{output_name}/magnify_{output_name}_markers.bam')
os.system(f'cnvpytor -root magnify_{output_name}/magnify.pytor -his {bin_size}')
os.system(f'cnvpytor -root magnify_{output_name}/magnify.pytor -partition {bin_size}')
os.system(f'cnvpytor -root magnify_{output_name}/magnify.pytor -call {bin_size} > magnify_{output_name}/calls.{bin_size}.tsv')
os.system(f'cnvpytor -root magnify_{output_name}/magnify.pytor -plot manhattan {bin_size} -o magnify_{output_name}/{output_name}.png')

os.system(f'samtools view -@ {threads} -H magnify_{output_name}/magnify_{output_name}_markers.bam > magnify_temp_{output_name}/magnify_{output_name}_markers_header.sam')

bed_file = region_bed(f"magnify_temp_{output_name}/magnify_{output_name}_markers_header.sam",marker_list)

os.system(f'samtools depth -@ {threads} -b {bed_file} magnify_{output_name}/magnify_{output_name}_markers.bam > magnify_temp_{output_name}/magnify_{output_name}_markers_coverage.cov')

chrHistograms(f'magnify_temp_{output_name}/magnify_{output_name}_markers_coverage.cov',marker_list)
os.system(f'mv fig_* magnify_{output_name}')
os.system(f'mv *.cov magnify_temp_{output_name}')

#KRAKEN2/BRACKEN system
if args.kraken[0].lower() == "t":
    os.system(f'samtools view -f 13 -@ {threads} magnify_{output_name}/magnify_{output_name}_markers.bam > magnify_temp_{output_name}/magnify_unaligned.sam')

    os.system(f'samtools fastq -@ {threads} magnify_temp_{output_name}/magnify_unaligned.sam -1 magnify_temp_{output_name}/magnify_unaligned_R1.fastq -2 magnify_temp_{output_name}/magnify_unaligned_R2.fastq')

    os.system(f'kraken2 --threads {threads} --db {database} --report magnify_{output_name}/classified_seqs_{output_name}.kreport --paired --classified-out magnify_{output_name}/classified_seqs_{output_name}#.fq magnify_temp_{output_name}/magnify_unaligned_R1.fastq magnify_temp_{output_name}/magnify_unaligned_R2.fastq > magnify_{output_name}/classified_{output_name}.kraken')

    os.system(f'bracken -d {database} -i magnify_{output_name}/classified_seqs_{output_name}.kreport -o magnify_{output_name}/classified_seqs_{output_name}.bracken -r 150')
