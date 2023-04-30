import os
import matplotlib.pyplot as plt

def region_bed(sam_header,chr_list):
    with open(f"{sam_header}_bed.bed","w") as bed:
        with open(f"{sam_header}") as file:
            for line in file:
                tags = line.split("\t")
                if tags[0] == '@SQ':
                    chr_name = tags[1].split(":")[1]
                    if chr_name in chr_list:
                        end = int(tags[2].split(":")[1])-1
                        bed.write(f"{chr_name}\t0\t{end}\n")
                else:
                    continue
    return f"{sam_header}_bed.bed"

def histogramData(coveragemap, chromosome, granularity=1):
    os.system(f"gawk '{{if ($1 ~ /({chromosome})\>/) print $0}};' {coveragemap} > {chromosome}_coverage.cov")
    with open(f"{chromosome}_coverage.cov") as file:
        bins, counts = [], []
        for line in file:
            features = line.split("\t")
            try:
                position, readnum = int(features[1].strip()), int(features[2].strip())
            except:
                continue
            try:
                if abs(bins[-1]-position) < granularity:
                    counts[-1] += readnum
                else:
                    bins.append(position)
                    counts.append(readnum)
            except:
                    bins.append(position)
                    counts.append(readnum)
    counts = [int(i)/granularity for i in counts]
    return bins,counts

def histogram(bins,counts,chromosome,granularity=50,rounding=1, significant=1):
    try:
        abnormal_bins, abnormal_counts = [], []

        mean = sum(counts)/len(counts)

        counts = [int(round((i/mean)*rounding,0)) for i in counts]

        mean_cnv = sum(counts)/len(counts)

        bins = [str(int(round(i/granularity,0))) for i in bins]
        
        std_dev = (sum([(i-mean_cnv)**2 for i in counts])/len(counts))**0.5
        for index in range(len(counts)):
            count = counts[index]
            if abs(count-mean_cnv) >= significant*std_dev:
                abnormal_bins.append(bins[index])
                abnormal_counts.append(counts[index])
                counts[index] = 0

        fig,ax = plt.subplots()
        ax.bar(bins,counts,color ='blue')
        ax.bar(abnormal_bins,abnormal_counts,color ='green')
        ax.axhline(y = mean_cnv, color = 'r', linestyle = '-')

        plt.xlabel(f"Base Coordinate (/{granularity})")
        plt.ylabel(f"Avg. Copy Number per Bin (x{rounding})")
        plt.xticks([])
        plt.title("Chromosome Coverage")
        fig.savefig(f'fig_{chromosome}.png')
        plt.close('all')
    except:
        pass

def chrHistograms(coveragemap, chrList):
    for i in range(len(chrList)):
        try:
            bins, counts = histogramData(coveragemap, chrList[i])
            histogram(bins, counts, chrList[i])
        except ZeroDivisionError:
            continue
