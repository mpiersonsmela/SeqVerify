import re
import os
import matplotlib.pyplot as plt
from math import exp, factorial, comb
from copy import deepcopy
from scipy.stats import poisson

supp_tags = ['SA','XA'] #sets the two possible optional alignments, chimeric and split respectively

def pathFinder(potential_path): #defines pathFinder, a function that checks if a variable is a filename or path
    if "/" in potential_path:
        path = potential_path
    else:
        path = os.getcwd()+"/"+potential_path
    return path

#def poisson(lam,x):
#    pdf = 0 #((lam**x)*(exp(-1*lam)))/(factorial(x))

#    to_sum = [(lam**j)/(factorial(j)) for j in range(0,x+1)]
#    cdf = (exp(-1*lam))*sum(to_sum)
#    return [pdf,cdf] 

#def binomial(n,p,k):
#    pdf = comb(n,k)*(p**k)*((1-p)**(n-k))
#    cdf = sum([comb(n,i)*(p**i)*((1-p)**(n-i)) for i in range(k)])
#    return [pdf,cdf]

class SamAlignment:
    def __init__(self, alignment): 
        #Defines all required parameters as specified in the SAM Manual to avoid dealing with slicing in future, as well as ALIGNMENT, a parameter containing the entire string, and OPTIONAL, a list containing all the optional fields
        
        #Whole alignment string
        self.ALIGNMENT = alignment

        #temp variable
        fields = alignment.split('\t')

        #SAM Headers
        self.QNAME = fields[0] 
        self.FLAG = fields[1] 
        self.RNAME = fields[2] 
        self.POS = fields[3]
        self.MAPQ = fields[4]
        self.CIGAR = fields[5]
        self.RNEXT = fields[6]
        self.PNEXT = fields[7]
        self.TLEN = fields[8]
        self.SEQ = fields[9]
        self.QUAL = fields[10]
        try:
            self.OPTIONAL = fields[11:]
        except IndexError:
            self.OPTIONAL = None #if optional fields do not exist, self.OPTIONAL returns None

    def __str__(self): #making an alignment into a string returns the original line from the sam file
        return self.ALIGNMENT

    def optionalTag(self,tag_code): #finds optional tag (if it exists, otherwise returns None) and separates it into [TAG,TYPE,VALUE] as determined in the SAM manual
        if self.OPTIONAL == None:
            return None
        for tag in self.OPTIONAL:
            if tag.startswith(tag_code):
                return tag.split(':')

    def supplementaryAlignments(self,tag_list=['SA','XA']): #finds supplementary alignments corresponding to the tags input (it expects a list). tags that follow this convention are SA (chimeric reads) and XA (split reads)
        supplementary_alignment_list = [] 
        for i in tag_list:
            optional_tag = self.optionalTag(i)
            if optional_tag == None:
                continue
            supplementary_alignment_string = optional_tag[2] #finds the SA tag (indicating 'other canonical alignments') and sets its value 
            supplementary_alignments = supplementary_alignment_string.split(';')
            supplementary_alignment_list = supplementary_alignment_list + [a.split(',') for a in supplementary_alignments if len(a) >= 2] #returns a list of lists of type [rname, pos, strand, CIGAR, mapQ, NM]
        return supplementary_alignment_list            

    def softClippingLen(self, cigar): #finds the length of the alignment until we hit softclipping in a string, if any, otherwise returns 0
        if "S" in cigar:
            cigar_split = re.split('(\d+)',cigar) #splits CIGAR string by regex matching groups of one or more digits
            cigar_letters = [cigar_split[i] for i in range(len(cigar_split)) if i % 2 == 0 and cigar_split[i] != ''] #takes the indicator letters of the CIGAR
            cigar_numbers = [cigar_split[i] for i in range(len(cigar_split)) if i % 2 == 1 and cigar_split[i] != ''] #takes the numbers after the letters
            length = 0
            if cigar_letters[0] == 'S': #if the first letter in the CIGAR denotes soft clipping, it means we do not need to adjust the read's length, so the function returns 0.
                return 0
            else:
                for operation_number in range(len(cigar_letters)):
                    letter = cigar_letters[operation_number]
                    number = int(cigar_numbers[operation_number])
                    if letter in ["M","I","S","=","X"]: #takes all the operations that change the length of the string, sums the length of their actions until we get to the clipping
                        if letter == "S":
                            return length
                        else:
                            length += number
                    else:
                        continue
        else:
            return 0

    def position(self):
        return [self.RNAME,int(self.POS)] #used to be just the int() statement

    def matePosition(self):
        return [self.RNEXT,int(self.PNEXT)+self.softClippingLen(self.optionalTag("MC")[-1])] #returns the position of a read's mate

    def supplementaryPosition(self,tag_list=['SA','XA']): #denotes where positions of supplementary alignments are exactly, given we know those are insertion sites
        positions = []
        supplementaries = self.supplementaryAlignments(tag_list) #finds all supplementary alignments
        for match in supplementaries:
            cigar = re.findall('(([0-9]+[A-Z])+)'," ".join(match[2:]))[0][0] #takes CIGAR string of the supplementary alignment
            position = match[1]
            soft_clipping_len = self.softClippingLen(cigar) #checks for softclipping
            if match[2] in ['-','+']: #match[1] and match[2] deal with conflicting ways in which position is reported, namely "-+NUM" in SA and ["-+","NUM"] in XA, so this makes sure we always get the correct position
                position = int(match[2] + str(position))

            if int(position) < 0: #adds the length of the string before any softclipping while respecting the orientation of the read
                position = abs(int(position) - soft_clipping_len) 
            else:
                position = int(position) + soft_clipping_len

            positions.append([match[0],-1*position]) #if a supplementary alignment is found, at the end we append its position multiplied by -1 to make sure it is not lost in granularity calculations with the non-supplementary alignments
        return positions

    def allPositions(self,tag_list=['SA','XA']): #returns position of the alignment, of its mate, and any supplementary alignments as a three-item list [[position],[mateposition],[supplementalaligmentposition]]
        required_positions = [self.position(),self.matePosition()]
        required_positions.extend(self.supplementaryPosition(tag_list))
        return required_positions

def group(samfile): #returns a dictionary containing all reference transgene chromosomes. each of those then contains another dictionary containing all the chromosomes that got mapped to, and each of those contains a dictionary containing the position of the mappings and how many times they were mapped to it.
    alignments = {}
    with open(samfile) as sam:
        for alignment in sam:
            if alignment.startswith("@"): #skips header lines
                continue
            alignment = SamAlignment(alignment) #uses our SamAlignment class for ease of use
            matches = alignment.allPositions(tag_list=supp_tags) #finds all supplementary alignments in the read
            ref_chromosome = matches[0][0] #takes name of ref chromosome
            mates_and_supplementaries = matches[1:] #takes positions of all matches and excludes the original read (we don't care where the read maps to on the original chromosome)
            try: #dictionary logic
                chromosome_dict = alignments[ref_chromosome] #initializes chromosome_dict to the reference chromosome entry, if it exists
                for position in mates_and_supplementaries:
                    aligned_chrm_name = position[0]
                    try:
                        try:
                            chromosome_dict[aligned_chrm_name][position[1]] += 1 
                        except:
                            chromosome_dict[aligned_chrm_name][position[1]] = 1 
                    except KeyError:
                        chromosome_dict[aligned_chrm_name] = {position[1]:1} 
            except KeyError:
                alignments[ref_chromosome] = {} #if the reference chromosome hasn't been seen yet, creates a new entry
                new_chromosomes = list(set([position[0] for position in mates_and_supplementaries])) #adds all the new chromosomes that have reads aligned to the reference chromosome in this alignment
                for chromosome in new_chromosomes: #initializes all of them to be empty
                    alignments[ref_chromosome][chromosome] = {}
                for position in mates_and_supplementaries: #if there is no dict for the chromosome the sequence aligns to, it creates one. this for loop is necessary as to not discard the information in that alignment
                    try:
                        alignments[ref_chromosome][position[0]][position[1]] += 1 
                    except:
                        alignments[ref_chromosome][position[0]][position[1]] = 1
    print("done with grouping")
    return alignments

def compress(alignment_dict, granularity=500): #compresses the alignments to the desired granularity
    print("reached compression")
    readout_dict = {} #initializes final dictionary as empty
    for reference_chromosome, alignments in alignment_dict.items():
        readout_dict[reference_chromosome] = {}
        for aligned_chromosome, alignment_data in alignments.items():
            readout_dict[reference_chromosome][aligned_chromosome] = {}
            for location, repetitions in alignment_data.items():
                saved_locations = readout_dict[reference_chromosome][aligned_chromosome].keys()
                distances = [abs(location - saved_location) for saved_location in saved_locations] #lines 148-154 cycle through the dictionary and copy its hierarchy, then finds the distances between a point and every other point for all points
                above_granularity = [distance >= granularity for distance in distances] #returns a list of true/false values whether a point is far enough that it needs to be in a different bin from the others
                if all(above_granularity): #if the point is far enough from all of them to be put into their bin, then...
                    try:
                        readout_dict[reference_chromosome][aligned_chromosome][location] += repetitions #...adds to a new bin or creates one
                    except KeyError:
                        readout_dict[reference_chromosome][aligned_chromosome][location] = repetitions
                else:
                    for possible_location, repetition in readout_dict[reference_chromosome][aligned_chromosome].items():
                        if abs(possible_location - location) < granularity: #if the point is at least close enough to one point to be put in its bin, finds which point it is and adds it to it.
                            readout_dict[reference_chromosome][aligned_chromosome][possible_location] += repetitions
    print("done with compression")
    return readout_dict

def filterAndScore(temp_folder,folder_insertion,bam_file,readout_dict,genome,read_depth):
    print("reached filtering")
    os.system(f'samtools depth {folder_insertion}/{bam_file} > {temp_folder}/total_coverage.cov')
    total_cov, length = 0,0
    with open(f"{temp_folder}/total_coverage.cov","r") as cov:
        for line in cov:
            depth = line.split("\t")[2]
            total_cov += int(depth)
            length += 1
    read_depth = (total_cov//length)
    print(f"calculated read depth as {read_depth} from {total_cov} total reads over length {length}")
    editable_readout = deepcopy(readout_dict)
    spurious_threshold = 0
    while poisson.cdf(spurious_threshold,read_depth) < 0.99999: #poisson(read_depth,spurious_threshold)[1]
        spurious_threshold += 1
    print(f"spurious threshold is {spurious_threshold}, proceeding to scoring...")
    with open(f"{temp_folder}/confidence.bed","w+") as bed:
        for read_chromosome, alignments in readout_dict.items():
            for alignment_chromosome, sites in alignments.items():
                for site, repetitions in sites.items():
                    repetitions = int(readout_dict[read_chromosome][alignment_chromosome][site])
                    bed.write(f"{alignment_chromosome}\t{abs(site)}\t{abs(int(site))+1}\n")
                    likelihood_insertion = poisson.cdf(repetitions,read_depth//2) #divides read length by 2 to get haploid depth #poisson(read_depth//2,repetitions)[1]
                    try:
                        ratio = likelihood_insertion/(likelihood_insertion+0.005) #TODO: Improve denominator
                        editable_readout[read_chromosome][alignment_chromosome][site] = round(ratio,4)
                    except ZeroDivisionError: 
                        editable_readout[read_chromosome][alignment_chromosome][site] = "inf"
    print("confidence score calculated, moving to coverage mapping...")
    os.system(f'samtools depth -b {temp_folder}/confidence.bed {folder_insertion}/{bam_file} > {temp_folder}/confidence.cov')
    with open(f"{temp_folder}/confidence.cov","r") as cov:
        to_delete = []
        for cov_site in cov:
            fields = cov_site.split("\t")
            if int(fields[2]) >= spurious_threshold:
                to_delete.append((fields[0],fields[1]))
    print("coverage mapping complete, moving to pruning high-coverage areas...")
    for read_chromosome, alignments in readout_dict.items():
        for alignment_chromosome, sites in alignments.items():
            for site, repetitions in sites.items():
                current_site = (alignment_chromosome,site)
                if current_site in to_delete:
                    del editable_readout[read_chromosome][alignment_chromosome][site]
    print("filtering complete")
    return [readout_dict,editable_readout]

def readout(folder,insertion_dict, original_dict, chr_filter, min_matches=1):
    print("reached printing")
    with open(f"{folder}/seqverify_readout.txt", "w") as file: #makes new readout file
        #file.write("Insertion Sites Found:\n") #boilerplate formatting
        file.write("chromosome,position,gene,chimeric,count,confidence\n")
        for read_chromosome, alignments in insertion_dict.items():
            print(read_chromosome, alignments)
            if alignments != {} and alignments != None:
                #file.write(read_chromosome+":\n") #inserts reference transgene first
                for align_chr, sites in alignments.items():
                    if sites != {}:
                        if align_chr in chr_filter:
                            continue #if the chromosome it aligns to is one we indicated we don't want (usually other transgenes, since due to repetitive DNA they clutter the readout), bins it and continues
                        else:
                            #file.write('\t'+align_chr+":\n") #adds the chromosome where insertion sites from the reference transgene were found
                            for site, repetitions in sites.items():
                                if repetitions != None:
                                    #if repetitions >= min_matches:
                                    match_num = original_dict[read_chromosome][align_chr][site]
                                    file.write(f"{align_chr},{str(site)[1:]},{read_chromosome},{'True' if str(site)[0] == '-' else 'False'},{str(match_num)},{str(repetitions)}\n")
                                    #if str(site)[0] == '-':
                                        #file.write('\t'+'\t'+str(site)[1:]+" (Chimeric/Split Read): "+str(repetitions)+" confidence score (matches)\n") #adds chimeric reads and their location
                                        #file.write(f"\t\t{str(site)[1:]} (Chimeric/Split Read): {str(repetitions)} confidence score ({str(match_num)} {'match' if int(match_num) == 1 else 'matches'}) \n") #adds chimeric reads and their location
                                    #else:
                                        #file.write('\t'+'\t'+str(site)+": "+str(repetitions)+" confidence score\n") #adds nonchimeric (potentially non-exaxt insertion sites) and their location
                                        #file.write(f"\t\t{str(site)}: {str(repetitions)} confidence score ({str(match_num)} {'match' if int(match_num) == 1 else 'matches'}) \n") #adds chimeric reads and their location
                                    #else:
                                    #    continue
    os.system(f'''awk 'NR<2 {{print $0;next}} {{print $0| "sort -t ',' -k3,3 -k1,1 -k2,2n "}}' {folder}/seqverify_readout.txt > {folder}/seqverify_readout.sorted.txt''')

def compare(vcf_1, vcf_2, min_quality, temp_folder, folder, stats, isec):
    for vcf in [vcf_1,vcf_2]:
        os.system(f"bgzip -f {vcf}")
        os.system(f"mv {vcf}.gz {temp_folder}")
        os.system(f"bcftools index {temp_folder}/{vcf}.gz")
    os.system(f"bcftools stats {temp_folder}/{vcf_1}.gz {temp_folder}/{vcf_2}.gz > {folder}/{stats}")
    id_dict = {'0':0,'1':0,'2':0}
    with open(f"{folder}/{stats}","r") as scores:
        for line in scores:
            if line.startswith("QUAL"):
                fields = line.split("\t")
                id, quality, freq = fields[1], fields[2], fields[3]
                if int(quality) >= min_quality:
                    id_dict[id] += int(freq)
    jaccard = str((id_dict['2'])/(id_dict['0']+id_dict['1']+id_dict['2']))
    with open(f"{folder}/{stats}","a") as scores:
        scores.write(f"The Jaccard similarity between {vcf_1} and {vcf_2} is {jaccard}")
    os.system(f"bcftools isec -p {temp_folder}/dir {temp_folder}/{vcf_1}.gz {temp_folder}/{vcf_2}.gz")
    os.system(f"mv {temp_folder}/dir/0001.vcf {folder}/{isec}")
