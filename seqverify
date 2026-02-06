#!/usr/bin/env python
 
import argparse
import os
import sys
import configparser
import shutil
import subprocess
import io
import warnings # Added to handle matplotlib warnings
from contextlib import redirect_stdout
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Import internal modules
from seqver_functions import *
from seqver_plots import *
from seqver_genomeupdate import *
from seqver_lofFinder import *
from seqver_gtf import *

# Import CNVpytor API
try:
    import cnvpytor
    from cnvpytor import Root, Genome, Viewer
except ImportError:
    print("Error: cnvpytor module not found. Please ensure it is installed.")
    sys.exit(1)

def run_cnvpytor_analysis(folder_cnv, temp_folder, aligned_bam, genome_fasta, output_name, bin_size, header_file):
    print("Initializing CNVpytor via Python API...")
    
    # 1. Clear default genomes to prevent resource check crash
    Genome.reference_genomes = {}
    
    # 2. Define Custom Genome from SAM Header
    chromosomes = OrderedDict()
    with open(header_file, "r") as f:
        for line in f:
            if line.startswith("@SQ"):
                parts = line.strip().split("\t")
                name = parts[1].split(":")[1]
                length = int(parts[2].split(":")[1])
                if name in ["chrM", "MT", "M"]:
                    chromosomes[name] = (length, "M")
                elif name in ["chrX", "X", "chrY", "Y"]:
                    chromosomes[name] = (length, "S")
                else:
                    chromosomes[name] = (length, "A")

    gc_file_path = os.path.join(temp_folder, f"{output_name}_gc.pytor")
    
    # Base configuration for custom genome
    custom_genome_conf = {
        "name": "custom",
        "species": "custom",
        "chromosomes": chromosomes
    }

    # 3. GC Correction Logic (Robust)
    gc_success = False
    if not os.path.exists(gc_file_path):
        print(f"Attempting to generate GC content file: {gc_file_path}")
        try:
            # Use a dedicated Root object to generate the GC file
            app_gc = Root(gc_file_path, create=True, max_cores=1)
            app_gc.gc(genome_fasta, make_gc_genome_file=True)
            gc_success = True
        except Exception as e:
            print(f"Warning: GC generation failed (proceeding without GC correction): {e}")
    else:
        gc_success = True

    if gc_success and os.path.exists(gc_file_path):
        print("GC file found/created. Enabling GC correction.")
        custom_genome_conf["gc_file"] = gc_file_path
    else:
        print("Skipping GC correction (file missing or generation failed).")

    # Update Global Config
    Genome.reference_genomes["custom"] = custom_genome_conf
    Genome.detected_genome = "custom"

    # 4. Initialize Root File
    pytor_file = os.path.join(folder_cnv, f"{output_name}.pytor")
    if os.path.exists(pytor_file):
        os.remove(pytor_file)
        
    app = Root(pytor_file, create=True)
    
    # 5. Run Steps
    print("Calculating Read Depth...")
    app.rd([aligned_bam], chroms=list(chromosomes.keys()))
    
    print(f"Calculating Histograms (bin={bin_size})...")
    app.calculate_histograms([bin_size])
    
    print("Partitioning...")
    app.partition([bin_size])
    
    print("Calling CNVs...")
    # 'call' method returns a dict of calls
    calls_dict = app.call([bin_size])
    
    # 6. Generate Call TSV
    call_file = os.path.join(folder_cnv, f"calls.{bin_size}.tsv")
    print(f"Writing calls to {call_file}...")
    with open(call_file, "w") as f:
        # Header matching CNVpytor output format
        f.write("type\tchrom\tstart\tend\tsize\tcnv\tp_val\tp_val_2\tp_val_3\tp_val_4\tQ0\tpN\tdG\n")
        if bin_size in calls_dict:
            for call in calls_dict[bin_size]:
                # Format call items to string and write tab-separated
                f.write("\t".join(map(str, call)) + "\n")
    
    # 7. Generate Plot
    print("Generating Manhattan Plot...")
    try:
        view = Viewer([pytor_file])
        
        # Configure params: disable GC correction if generation failed
        view.params.update({
            "bin_size": bin_size,
            "chrom": list(chromosomes.keys()),
            "rd_use_gc_corr": gc_success, 
            "rd_manhattan_call": True 
        })
        
        plt.figure(figsize=(12, 6))
        
        # [FIX] Suppress matplotlib warnings caused by CNVpytor re-initializing figures
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, message="Ignoring specified arguments")
            warnings.filterwarnings("ignore", category=UserWarning, message="FigureCanvasAgg is non-interactive")
            view.manhattan() 
            
        plot_path = os.path.join(folder_cnv, f"{output_name}.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()
        print(f"Plot saved to {plot_path}")
    except Exception as e:
        print(f"Error during plotting: {e}")

# --- MAIN SCRIPT START ---

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=str, required=False, default="output") 
    parser.add_argument('--reads_1', type=str, required=False) 
    parser.add_argument('--reads_2', type=str, required=False) 
    parser.add_argument('--untargeted', type=str, required=False, nargs='+') 
    parser.add_argument('--targeted', type=str, required=False) 
    parser.add_argument('--keepgoing',action='store_const',const=True,required=False,default=False) 
    parser.add_argument('--keep_temp', action='store_const',const=True, required=False, default=False) 
    parser.add_argument('--use_mem2', action='store_const',const=True, required=False,default=False) 

    parser.add_argument('--genome', type=str, required=False, default=f"{os.getcwd()}/seqverify_defaults/chm13v2.0.fa") 
    parser.add_argument('--gtf', type=str, required=False, default=f"{os.getcwd()}/seqverify_defaults/chm13v2.0_RefSeq_Liftoff_v5.1.gff3") 

    parser.add_argument('--kraken', action='store_const',const=True, required=False,default=False) 
    parser.add_argument('--database',type=str,required=False, default=f"{os.getcwd()}/seqverify_defaults/seqverify_database") 

    parser.add_argument('--granularity',type=int, required=False, default=500) 
    parser.add_argument('--threads',type=int,required=False, default=1) 
    parser.add_argument('--max_mem', type=str, required=False, default='16G') 
    parser.add_argument('--min_matches', type=int, required=False,default=1) 
    parser.add_argument('--start', type=str, required=False, default="all") 
    parser.add_argument('--mitochondrial',action='store_const',const=True,required=False,default=False) 
    parser.add_argument('--stringency',type=float, required=False, default=0.005) 
    parser.add_argument('--spurious_filtering_threshold', type=float, required=False, default=0.00001)

    parser.add_argument('--manual_plots',action='store_const',const=True,required=False,default=False) 
    parser.add_argument('--bin_size',type=int, required=False, default=100000) 

    parser.add_argument('--variant_calling',nargs=2,type=str,required=False) 
    parser.add_argument('--variant_intensity',type=str,required=False,default="MODERATE") 
    parser.add_argument('--variant_window_size',type=int,required=False,default=10000) 

    parser.add_argument('--download_defaults',action='store_const',const=True,required=False,default=False) 
    parser.add_argument('--similarity',nargs=2,type=str,required=False) 
    parser.add_argument('--min_quality',type=int,required=False,default=3) 
    parser.add_argument('--config',type=str,required=False,default=None)

    args = parser.parse_args()

    if args.config is not None:
        config = configparser.ConfigParser(allow_no_value=True)
        config.read(args.config)
        args.output = config["MAIN"]["output"]
        args.reads_1 = config["MAIN"]["reads_1"]
        args.reads_2 = config["MAIN"]["reads_2"]
        args.untargeted = config["MAIN"]["untargeted"].split(" ")
        args.targeted = config["MAIN"]["targeted"]
        if not args.targeted.endswith(".txt"): args.targeted = None 
        args.keepgoing = bool(config["MAIN"]["keepgoing"] == "True")
        args.keep_temp = bool(config["MAIN"]["keep_temp"] == "True")
        args.use_mem2 = bool(config["MAIN"]["use_mem2"] == "True")
        args.genome = config["GENOME"]["genome"]
        args.gtf = config["GENOME"]["gtf"]
        args.kraken = bool(config["KRAKEN2"]["kraken"] == "True")
        args.database = config["KRAKEN2"]["database"]
        args.granularity = int(config["INSERTION"]["granularity"])
        args.threads = int(config["INSERTION"]["threads"])
        args.max_mem = config["INSERTION"]["max_mem"]
        args.min_matches = int(config["INSERTION"]["min_matches"])
        args.start = config["INSERTION"]["start"]
        args.mitochondrial = bool(config["INSERTION"]["mitochondrial"] == "True")
        args.stringency = float(config["INSERTION"]["stringency"])
        args.spurious_filtering_threshold = float(config["INSERTION"]["spurious_filtering_threshold"])
        args.manual_plots = bool(config["CNV"]["manual_plots"] == "True")
        args.bin_size = int(config["CNV"]["bin_size"])
        args.variant_calling = eval(config["VARIANT"]["variant_calling"])
        args.variant_intensity = config["VARIANT"]["variant_intensity"]
        args.variant_window_size = int(config["VARIANT"]["variant_window_size"])
        args.min_quality = int(config["OTHER"]["min_quality"])

    print("Running SeqVerify with arguments:")
    for arg, value in vars(args).items():
        print(f"{arg}: {value}")

    if args.download_defaults: 
        os.system('mkdir seqverify_defaults')
        os.system('curl -OJX GET "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz" --output-dir seqverify_defaults')
        os.system('gunzip seqverify_defaults/chm13v2.0.fa.gz')
        os.system('curl -OJX GET "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz" --output-dir seqverify_defaults')
        os.system('gunzip seqverify_defaults/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz')
        os.system('curl -OJX GET "ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" --output-dir seqverify_defaults')
        os.system('gunzip seqverify_defaults/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
        os.system('mkdir seqverify_defaults/seqverify_database')
        os.system('curl -OJX GET "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08gb_20230605.tar.gz" --output-dir seqverify_defaults')
        os.system('tar -xzf seqverify_defaults/k2_pluspf_08gb_20230605.tar.gz -C seqverify_defaults/seqverify_database')
        os.system('curl -OJX GET "https://raw.githubusercontent.com/pcingola/SnpEff/master/config/snpEff.config" --output-dir seqverify_defaults')
        os.system('curl -OJX GET "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" --output-dir seqverify_defaults')
        os.system('curl -OJX GET "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi" --output-dir seqverify_defaults')
        os.system('gunzip seqverify_defaults/clinvar.vcf.gz')
        sys.exit(0)

    output_name = args.output
    folder = f"{output_name}_seqverify"
    temp_folder =f"{output_name}_seqverify_temp"
    marker_list_file = f"seqverify_{output_name}_transgene_list.fa"
    genome_with_markers = f"seqverify_{output_name}_genome_with_markers.fa"
    aligned_with_markers = f"seqverify_{output_name}_markers.sam"
    aligned_diff_chr = f"seqverify_{output_name}_markers_diff_chr.sam"
    aligned_diff_chr_bam = f"seqverify_{output_name}_markers_diff_chr.bam"
    header = f"seqverify_{output_name}_markers_header.sam"
    overall_coverage = f"seqverify_{output_name}_markers_coverage.cov"
    pytor = f"{output_name}.pytor"
    unaligned = f"seqverify_{output_name}_unaligned.sam"
    unaligned_R1 = f"seqverify_{output_name}_unaligned_R1.fastq"
    unaligned_R2 = f"seqverify_{output_name}_unaligned_R2.fastq"
    variant_genome_aligned = f"seqverify_{output_name}_variant.sam"
    variant_genome_aligned_bam = f"seqverify_{output_name}_variant.bam"
    variant_vcf = f"seqverify_{output_name}.vcf"
    variant_vcf_ann = f"seqverify_{output_name}.ann.vcf"
    variant_lof = f"seqverify_{output_name}_variants.tsv"
    bed_file = f"{output_name}.bed"
    vcf_similarity = f"seqverify_{output_name}_bcfstats.txt"
    vcf_isec = f"seqverify_{output_name}_unique.vcf"
    new_gtf_unsorted = f"seqverify_{output_name}_unsorted.gtf"
    new_gtf_sorted = f"seqverify_{output_name}.gtf"

    folder_insertion = f"{folder}/insertion"
    folder_cnv = f"{folder}/copy_number"
    folder_kraken = f"{folder}/kraken"
    folder_snp =f"{folder}/variant_calling"

    if args.similarity is not None:
        compare(args.similarity[0],args.similarity[1],args.min_quality,temp_folder,folder,vcf_similarity, vcf_isec)
        sys.exit(0)

    reads_1_path = pathFinder(args.reads_1)
    reads_2_path = pathFinder(args.reads_2)
    genome_source_path = pathFinder(args.genome)
    marker_sources_path = None
    exact_path = None
    database = None
    gtf_file = None

    granularity_threshold = args.granularity

    if args.untargeted is not None and args.untargeted != ['']:
        print("Loading " + str(len(args.untargeted)) + " untargeted insertions:" + str(args.untargeted))
        marker_sources_path = [pathFinder(i) for i in args.untargeted]
    else:
        print("No untargeted insertions loaded")

    if args.targeted is not None:
        print("Loading targeted insertions:" + args.targeted)
        exact_path = pathFinder(args.targeted)
    else:
        print("No targeted insertions loaded")

    if args.database is not None:
        database = pathFinder(args.database)
    if args.gtf is not None:
        gtf_file = pathFinder(args.gtf)
        
    os.system(f'mkdir -p {folder}')
    os.system(f'mkdir -p {temp_folder}')    

    if exact_path is not None:
        commands = commandHandler(genome_source_path,exact_path,temp_folder,genome_with_markers,return_commands_only = False)

    min_matches = args.min_matches
    bin_size = args.bin_size
    threads = args.threads

    if args.bin_size != 100000:
        if args.bin_size % 100 == 0:
            bin_size = args.bin_size
        else:
            raise ValueError("Custom bin size not divisible by 100")

    if args.max_mem[-1] == "M":
        max_mem = int(args.max_mem[:-1])*1000000
    elif args.max_mem[-1] == "G":
        max_mem = int(args.max_mem[:-1])*1000000000
    elif type(args.max_mem) == 'int':
        max_mem = int(args.max_mem)

    base_block = max_mem / 8
    os.environ['_JAVA_OPTIONS'] = '-Xmx'+args.max_mem

    if args.start == "beginning" or args.start == "all":
        os.system(f'mkdir -p {folder_insertion}')
        print("started gtf")
        if args.gtf is not None and exact_path is not None:
            gtfEdit(gtf_file, commands, temp_folder, folder_insertion, new_gtf_unsorted, new_gtf_sorted)
        
        if os.path.exists(f'{temp_folder}/{marker_list_file}'):
            os.remove(f'{temp_folder}/{marker_list_file}')
        os.system(f'touch {temp_folder}/{marker_list_file}')

        if exact_path is None:
            os.system(f'cp {genome_source_path} {temp_folder}/{genome_with_markers}')

        if marker_sources_path is not None:
            for file in marker_sources_path: 
                print(f'Reading marker file: {file}')
                os.system(f'cat {file} >> {temp_folder}/{marker_list_file}') 

            print(f'Appending markers to genome...')
            os.system(f'echo "" >> {temp_folder}/{genome_with_markers}')
            os.system(f'cat {temp_folder}/{marker_list_file} >> {temp_folder}/{genome_with_markers}')

        if args.use_mem2:
            os.system(f'bwa-mem2 index {temp_folder}/{genome_with_markers}')
        else:
            os.system(f'bwa index -b {base_block} {temp_folder}/{genome_with_markers}')
        os.system(f'cp {temp_folder}/{genome_with_markers} {folder_insertion}/{genome_with_markers}')
        
        if args.keepgoing:
            print("Finished beginning, continuing to align")
            args.start = "align"

    if args.start == "align" or args.start == "all":
        if args.use_mem2:
            os.system(f'bwa-mem2 mem -M -t {threads} {temp_folder}/{genome_with_markers} {reads_1_path} {reads_2_path} > {temp_folder}/{aligned_with_markers}')
        else:
            os.system(f'bwa mem -M -t {threads} {temp_folder}/{genome_with_markers} {reads_1_path} {reads_2_path} > {temp_folder}/{aligned_with_markers}')
        if args.keepgoing:
            print("Finished align, continuing to markers")
            args.start = "markers"
        
    if args.start == "markers" or args.start == "all":
        os.system(f'mkdir -p {folder_insertion}')
        
        has_markers = False
        if os.path.exists(f'{temp_folder}/{marker_list_file}'):
            with open(f'{temp_folder}/{marker_list_file}') as f:
                marker_list = [i[1:].strip() for i in f if i.startswith('>')]
                if marker_list:
                    has_markers = True
                    markers = "|".join(marker_list)
        
        if has_markers:
            if args.mitochondrial:
                markers += "|chrM"
                os.system(f"samtools view -h -@ {threads} {temp_folder}/{aligned_with_markers} | awk '{{if ((($3 ~ /^({markers})/)&&($7 ~ /chr([0-9]+|[XY])/))||($1 ~ /^@/)) print $0}};' > {temp_folder}/{aligned_diff_chr}")
            else:
                os.system(f"samtools view -h -@ {threads} {temp_folder}/{aligned_with_markers} | awk '{{if ((($3 ~ /^({markers})/)&&($7 ~ /chr([0-9]+|[XYM])/))||($1 ~ /^@/)) print $0}};' > {temp_folder}/{aligned_diff_chr}")
        else:
            os.system(f'samtools view -H {temp_folder}/{aligned_with_markers} > {temp_folder}/{aligned_diff_chr}')

        os.system(f"samtools sort -@ {threads} {temp_folder}/{aligned_with_markers} > {folder_insertion}/{aligned_diff_chr_bam}")
        os.system(f"samtools index -@ {threads} {folder_insertion}/{aligned_diff_chr_bam}")
        os.system(f'samtools view -@ {threads} -H {folder_insertion}/{aligned_diff_chr_bam} > {temp_folder}/{header}')
        
        if args.keepgoing:
            print("Finished markers, continuing to readout")
            args.start = "readout"
            
    if args.start == "readout" or args.start == "all":
        marker_list = []
        if os.path.exists(f'{temp_folder}/{marker_list_file}'):
            with open(f'{temp_folder}/{marker_list_file}') as f:
                marker_list = [i[1:].strip() for i in f if i.startswith('>')]
        
        data = group(f'{temp_folder}/{aligned_diff_chr}')
        insertions = compress(data,granularity_threshold)

        if insertions and len(insertions) > 0:
            filteredInsertions = filterAndScore(temp_folder,folder_insertion,aligned_diff_chr_bam,insertions,args.spurious_filtering_threshold,args.stringency)
            readout(folder_insertion,filteredInsertions[1],filteredInsertions[0],marker_list,min_matches)
        else:
            print("No insertions found, skipping filtering/scoring.")
            readout(folder_insertion, {}, {}, marker_list, min_matches)
        
        if args.keepgoing:
            print("Finished markers, continuing to CNV")
            args.start = "cnv"

    if args.start == "cnv" or args.start == "all":
        os.system(f'mkdir -p {folder_cnv}')
        
        run_cnvpytor_analysis(
            folder_cnv=folder_cnv,
            temp_folder=temp_folder,
            aligned_bam=f"{folder_insertion}/{aligned_diff_chr_bam}",
            genome_fasta=f"{temp_folder}/{genome_with_markers}",
            output_name=output_name,
            bin_size=bin_size,
            header_file=f"{temp_folder}/{header}"
        )

        if args.keepgoing:
            print("Finished CNVpytor, continuing to plots")
            args.start = "plots"

    if args.start == "plots" or args.start == "all":
        os.system(f'mkdir -p {folder_cnv}')
        
        if 'marker_list' not in locals():
            marker_list = []
            if os.path.exists(f'{temp_folder}/{marker_list_file}'):
                with open(f'{temp_folder}/{marker_list_file}') as f:
                    marker_list = [i[1:].strip() for i in f if i.startswith('>')]

        try:
            region_bed(temp_folder,header,commands,marker_list,bed_file)
        except NameError:
            if exact_path is not None:
                commands = commandHandler(genome_source_path,exact_path,temp_folder,genome_with_markers, return_commands_only = True)
            else:
                commands = None
            
            markers = "|".join(marker_list)
            region_bed(temp_folder,header,commands,marker_list,bed_file)
            
        if not args.manual_plots:
            print("Plotting transgenes with igv-reports")
            if args.gtf is not None:
                if exact_path is not None:
                    gtf = new_gtf_sorted 
                else:
                    if os.path.exists(args.gtf):
                        gtf_basename = os.path.basename(args.gtf)
                        shutil.copy(args.gtf, f"{folder_insertion}/{gtf_basename}")
                        gtf = gtf_basename
                    else:
                        gtf = args.gtf
            igvScreenshot_new(temp_folder,folder_insertion,f'{folder_insertion}/{aligned_diff_chr_bam}',f'{folder_insertion}/{genome_with_markers}',bed_file, gtf_file=gtf)
        
        # [FIX] Force manual plotting of ALL chromosomes (including main genome)
        # This provides a clean "second opinion" plot using raw coverage data
        # regardless of whether the user asked for it or not, which is helpful for verification.
        print("Generating Full Genome Coverage Plots (Manual Fallback)...")
        
        # 1. Get all chromosome names from the BAM header
        all_chroms = []
        with open(f"{temp_folder}/{header}", "r") as f:
            for line in f:
                if line.startswith("@SQ"):
                    # Extract "SN:chr1" -> "chr1"
                    all_chroms.append(line.split("\t")[1].split(":")[1])
        
        # 2. Run coverage calculation for everything
        # We reuse the 'overall_coverage.cov' file if it was generated, or generate it now
        cov_file = f'{temp_folder}/{overall_coverage}'
        # Using a simple 0-based BED for the whole genome is cleaner than region_bed for this purpose
        # But we can just run depth on the whole BAM without regions
        os.system(f'samtools depth {folder_insertion}/{aligned_diff_chr_bam} > {cov_file}')
        
        # 3. Generate histograms
        # Note: chrHistograms was originally designed for markers, but works for any chrom
        # in the list.
        chrHistograms(cov_file, all_chroms)
        
        # 4. Move plots
        os.system(f'mv fig_* {folder_insertion}')
        os.system(f'mv *.cov {temp_folder}')
            
        if args.keepgoing:
            print("Finished plots, continuing to Kraken (if enabled)")
            args.start = "kraken"

    if args.start == "kraken" or args.start == "all":
        os.system(f'mkdir -p {folder_kraken}')

        if args.kraken:
            os.system(f'samtools view -f 13 -@ {threads} {folder_insertion}/{aligned_diff_chr_bam}> {temp_folder}/{unaligned}')
            os.system(f'samtools fastq -@ {threads} {temp_folder}/{unaligned} -1 {temp_folder}/{unaligned_R1} -2 {temp_folder}/{unaligned_R2}')
            os.system(f'kraken2 --threads {threads} --db {database} --report {folder_kraken}/classified_seqs_{output_name}.kreport --paired --classified-out {folder_kraken}/classified_seqs_{output_name}#.fq {temp_folder}/{unaligned_R1} {temp_folder}/{unaligned_R2} > {folder_kraken}/classified_{output_name}.kraken')
            os.system(f'bracken -d {database} -i {folder_kraken}/classified_seqs_{output_name}.kreport -o {folder_kraken}/classified_seqs_{output_name}.bracken -r 150')

        if args.keepgoing:
            print("Finished kraken, continuing to variant calling (if enabled)")
            args.start = "variant"

    if args.start == "variant" or args.start == "all":
        os.system(f'mkdir -p {folder_snp}')
        if args.variant_calling is not None and args.variant_calling != []:
            variant_genome_source_path = pathFinder(args.variant_calling[0]) 
            clinvardb_source_path = pathFinder(args.variant_calling[1])
            
            if args.use_mem2:
                print("Creating BWA-MEM2 index")
                os.system(f"bwa-mem2 index {variant_genome_source_path}")
            else:
                if os.path.isfile(variant_genome_source_path+".bwt") and os.path.isfile(variant_genome_source_path+".pac") and os.path.isfile(variant_genome_source_path+".sa") and os.path.isfile(variant_genome_source_path+".ann") and os.path.isfile(variant_genome_source_path+".amb"):
                    print("Detected that BWA indexing was already performed. Skipping BWA indexing.")
                else:
                    print("Existing BWA index not detected. Performing indexing.")
                    os.system(f"bwa index {variant_genome_source_path}")
            
            if os.path.isfile(variant_genome_source_path+".fai"):
                print("Detected that samtools faidx was already performed. Skipping samtools faidx.")
            else:
                print("Existing FASTA index not detected. Performing indexing.")
                os.system(f"samtools faidx {variant_genome_source_path}")
            
            if args.use_mem2:
                bwa_cmd = f"bwa-mem2 mem -M -t {threads} {variant_genome_source_path} {reads_1_path} {reads_2_path} > {temp_folder}/{variant_genome_aligned}"
            else:
                bwa_cmd = f"bwa mem -M -t {threads} {variant_genome_source_path} {reads_1_path} {reads_2_path} > {temp_folder}/{variant_genome_aligned}"
            
            exitcode = os.system(bwa_cmd)
            if exitcode!=0: sys.exit("BWA failed")
            
            os.system(f"samtools sort -@ {threads} {temp_folder}/{variant_genome_aligned} > {temp_folder}/{variant_genome_aligned_bam}") 
            os.system(f"samtools index -@ {threads} {temp_folder}/{variant_genome_aligned_bam}")
            os.system(f"bcftools mpileup -Ou -f {variant_genome_source_path} {temp_folder}/{variant_genome_aligned_bam} | bcftools call -mv -Ov -o {temp_folder}/{variant_vcf}")
            
            if args.keepgoing:
                args.start = "snp_filtering"
            
    if (args.start == "snp_filtering" or args.start == "all") and args.variant_calling is not None:
            variant_genome_source_path = pathFinder(args.variant_calling[0])
            clinvardb_source_path = pathFinder(args.variant_calling[1])
            
            os.system(f"snpEff -dataDir seqverify_defaults -v GRCh38.105 {temp_folder}/{variant_vcf} > {temp_folder}/{variant_vcf_ann}")
            os.system(f"SnpSift annotate -v {clinvardb_source_path} {temp_folder}/{variant_vcf_ann} > {folder_snp}/{variant_vcf_ann}")
            
            if exact_path is None:
                mutation_logger(folder_snp,variant_lof,variant_vcf_ann,args.min_quality,args.variant_intensity, args.variant_window_size)
            else:
                mutation_logger(folder_snp,variant_lof,variant_vcf_ann,args.min_quality,args.variant_intensity, args.variant_window_size, exact_path)

    if not args.keep_temp:
        os.system(f'rm -r {temp_folder}')

if __name__ == "__main__":
    main()