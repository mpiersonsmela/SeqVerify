#Defines functions necessary for the plotting of transgene copy numbers, either through IGV or matplotlib.
import os
import matplotlib.pyplot as plt

def region_bed(temp_folder,sam_header,chr_list): #Generates .bed file necessary for running the other functions
    with open(f"{temp_folder}/{sam_header}_bed.bed","w") as bed: #Creates a new bed file
        with open(f"{sam_header}") as file:
            for line in file:
                tags = line.split("\t") #Parses the file line-by-line, finding each tab-separated tag.
                if tags[0] == '@SQ': #Collects all sequence/chromosome lines in the header
                    chr_name = tags[1].split(":")[1] #Takes the chromosome/transgene's name
                    if chr_name in chr_list: #Excludes all the names we don't want (i.e. usually all the human non-transgene chromosomes)
                        end = int(tags[2].split(":")[1])-1 #Calculates 0-indexed final coordinate of transgene (SAM files are 1-indexed, .bed files are )
                        bed.write(f"{chr_name}\t0\t{end}\n") #Writes the coordinates to the file
                else:
                    continue
    return f"{sam_header}_bed.bed"

def histogramData(coveragemap, chromosome, granularity=1): #Collects data to make a single histogram if IGV is not used
    os.system(f"gawk '{{if ($1 ~ /({chromosome})\>/) print $0}};' {coveragemap} > {chromosome}_coverage.cov") #gawks the output of samtools depth -b for the reads relevant to the chromosome
    with open(f"{chromosome}_coverage.cov") as file:
        bins, counts = [], [] #initializes an empty list of bins and an empty list of read numbers per bin 
        for line in file:
            features = line.split("\t")
            try:
                position, readnum = int(features[1].strip()), int(features[2].strip()) #takes the position and number of reads at that position for each line in the file
            except:
                continue
            try:
                if abs(bins[-1]-position) < granularity: #if the bin that was just analyzed falls within a previous bin, add the number of reads to the previous bin
                    counts[-1] += readnum
                else: #if not, make a new bin and a new corresponding number of reads
                    bins.append(position)
                    counts.append(readnum)
            except:
                    bins.append(position)
                    counts.append(readnum)
    counts = [int(i)/granularity for i in counts] #averages the readnums over the length of the bins as to indicate average read number
    return bins,counts

def histogram(bins,counts,chromosome,granularity=50,rounding=1, significant=1): #generates the plot if using matplotlib
    try:
        abnormal_bins, abnormal_counts = [], [] #initializes bins and read counts for any potentially abnormal bins

        mean = sum(counts)/len(counts) #calculates the average number of reads in the transgene overall

        counts = [int(round((i/mean)*rounding,0)) for i in counts] #calculates the mean copy number for all the bins 

        mean_cnv = sum(counts)/len(counts) #calculates the mean copy number overall

        bins = [str(int(round(i/granularity,0))) for i in bins] #calculates the new bins
        
        std_dev = (sum([(i-mean_cnv)**2 for i in counts])/len(counts))**0.5 #calculates the overall standard deviation in the copy numbers
        for index in range(len(counts)):
            count = counts[index]
            if abs(count-mean_cnv) >= significant*std_dev: #if the difference between a copy number and the mean copy number is larger than the threshold coefficient multiplied by the standard deviation, flags the bin as abnormal
                abnormal_bins.append(bins[index])
                abnormal_counts.append(counts[index])
                counts[index] = 0

        fig,ax = plt.subplots() #matplotlib logic
        ax.bar(bins,counts,color ='blue') #plots all normal bins in blue
        ax.bar(abnormal_bins,abnormal_counts,color ='green') #plots all abnormal/significant bins in green
        ax.axhline(y = mean_cnv, color = 'r', linestyle = '-') #plots a red average line

        #Defines labels, saves the image, closes the plot
        plt.xlabel(f"Base Coordinate (/{granularity})") 
        plt.ylabel(f"Avg. Copy Number per Bin (x{rounding})")
        plt.xticks([])
        plt.title("Chromosome Coverage")
        fig.savefig(f'fig_{chromosome}.png')
        plt.close('all')
    except:
        pass

def chrHistograms(coveragemap, chrList): #Generates histograms for every chromosome in chrList by calling the above two functions each time.
    for i in range(len(chrList)):
        try:
            bins, counts = histogramData(coveragemap, chrList[i])
            histogram(bins, counts, chrList[i])
        except ZeroDivisionError:
            continue

def igvScreenshot(temp_folder,folder,alignments,genome,bed_file): #If using IGV, deals with IGV logic
    with open(f"{temp_folder}/seqverify_igv.bat","w+") as file: #Autogenerates an IGV-compatible bat file we will use later to screenshot the relevant parts
        file.write("new\n") #boilerplate code
        file.write(f"snapshotDirectory {folder}\n") #sets the folder for the screenshots
        file.write(f"load {alignments}\n") #loads the bam file
        file.write(f"genome {genome}\n") #loads the genome
        file.write(f"maxPanelHeight 500\n") #boilerplate code for adjustment of the screen size
        with open(f"{temp_folder}/{bed_file}","w+") as bed: #writes instructions to take a screenshot of every transgene in the bed file
            for line in bed:
                name,begin,end = line.split("\t")
                file.write(f"goto {name}:{begin}-{end}") #makes IGV go to the entire transgene
                file.write(f"snapshot fig_{name}.png") #takes screenshot and saves it
        file.write("exit") #boilerplate
    os.system("xvfb-run --auto-servernum --server-num=1 igv -b seqverify_igv.bat") #runs XVFB, a headerless server emulator, to run IGV automatically without the need for a GUI.
