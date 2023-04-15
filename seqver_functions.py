import regex as re
import os

supp_tags = ['SA','XA']
chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
#grep A00738:415:HLNKWDSX3:2:2608:7193:18490 chm13v2.0_realigned_transgenes.sam
#python magnify.py --reads_source S04_sorted_marked.bam --genome_source chm13v2.0.fa --marker_sources T2A-mGreenLantern.fa T2A-tdTomato.fa transposons_in_S04.fa unwanted_plasmids.fa --granularity 500 --threads 15
def pathFinder(potential_path):
    if "/" in potential_path:
        path = potential_path
    else:
        path = os.getcwd()+"/"+potential_path
    return path

class SamAlignment:
    def __init__(self, alignment): #defines all required parameters as specified in the SAM Manual to avoid dealing with slicing in future, as well as ALIGNMENT, a parameter containing the entire string, and OPTIONAL, a list containing all the optional fields

        self.ALIGNMENT = alignment

        fields = alignment.split('\t')

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

    def softClippingLen(self, cigar):
        if "S" in cigar:
            cigar_split = re.split('(\d+)',cigar)
            cigar_letters = [cigar_split[i] for i in range(len(cigar_split)) if i % 2 == 0 and cigar_split[i] != ''] 
            cigar_numbers = [cigar_split[i] for i in range(len(cigar_split)) if i % 2 == 1 and cigar_split[i] != '']
            length = 0
            if cigar_letters[0] == 'S':
                return 0
            else:
                for operation_number in range(len(cigar_letters)):
                    letter = cigar_letters[operation_number]
                    number = int(cigar_numbers[operation_number])
                    if letter in ["M","I","S","=","X"]:
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
        return [self.RNEXT,int(self.PNEXT)+self.softClippingLen(self.optionalTag("MC")[-1])]

    def supplementaryPosition(self,tag_list=['SA','XA']): #you minus one'd the coords here to designate them
        positions = []
        supplementaries = self.supplementaryAlignments(tag_list)
        for match in supplementaries:
            cigar = re.findall('(([0-9]+[A-Z])+)'," ".join(match[2:]))[0][0]
            position = match[1]
            soft_clipping_len = self.softClippingLen(cigar)
            if match[2] in ['-','+']:
                position = int(match[2] + str(position))

            if int(position) < 0:
                position = abs(int(position) - soft_clipping_len) #used to be a minus
            else:
                position = int(position) + soft_clipping_len

            positions.append([match[0],-1*position])
        return positions

    def allPositions(self,tag_list=['SA','XA']): #returns position of the alignment, of its mate, and any supplementary alignments as a three-item list [[position],[mateposition],[supplementalaligmentposition]]
        required_positions = [self.position(),self.matePosition()]
        required_positions.extend(self.supplementaryPosition(tag_list))
        return required_positions

def group(samfile): #returns a dictionary containing all reference transgene chromosomes. each of those then contains another dictionary containing all the chromosomes that got mapped to, and each of those contains a dictionary containing the position of the mappings and how many times they were mapped to it.
    alignments = {}
    with open(samfile) as sam:
        for alignment in sam:
            if alignment.startswith("@"):
                continue
            alignment = SamAlignment(alignment)
            matches = alignment.allPositions(tag_list=supp_tags)
            ref_chromosome = matches[0][0] #takes name of ref chromosome
            mates_and_supplementaries = matches[1:] #takes positions of all matches and excludes the original read (we don't care where the read maps to on the original chromosome)
            try:
                chromosome_dict = alignments[ref_chromosome]
                for position in mates_and_supplementaries:
                    aligned_chrm_name = position[0]
                    try:
                        try:
                            chromosome_dict[aligned_chrm_name][position[1]] += 1 #there used to be an abs here
                        except:
                            chromosome_dict[aligned_chrm_name][position[1]] = 1 #abs here
                    except KeyError:
                        chromosome_dict[aligned_chrm_name] = {position[1]:1}  #abs here
            except KeyError:
                alignments[ref_chromosome] = {}
                new_chromosomes = list(set([position[0] for position in mates_and_supplementaries]))
                for chromosome in new_chromosomes:
                    alignments[ref_chromosome][chromosome] = {}
                for position in mates_and_supplementaries: #if there is no dict for the chromosome the sequence aligns to, it creates one. this for loop is necessary as to not discard the information in that alignment
                    try:
                        alignments[ref_chromosome][position[0]][position[1]] += 1
                    except:
                        alignments[ref_chromosome][position[0]][position[1]] = 1
    return alignments

def compress(alignment_dict, granularity=500):
    readout_dict = {}
    for reference_chromosome, alignments in alignment_dict.items():
        readout_dict[reference_chromosome] = {}
        for aligned_chromosome, alignment_data in alignments.items():
            readout_dict[reference_chromosome][aligned_chromosome] = {}
            for location, repetitions in alignment_data.items():
                saved_locations = readout_dict[reference_chromosome][aligned_chromosome].keys()
                distances = [abs(location - saved_location) for saved_location in saved_locations]
                above_granularity = [distance >= granularity for distance in distances]
                if all(above_granularity):
                    try:
                        readout_dict[reference_chromosome][aligned_chromosome][location] += repetitions
                    except KeyError:
                        readout_dict[reference_chromosome][aligned_chromosome][location] = repetitions
                else:
                    for possible_location, repetition in readout_dict[reference_chromosome][aligned_chromosome].items():
                        if abs(possible_location - location) < granularity:
                            readout_dict[reference_chromosome][aligned_chromosome][possible_location] += repetitions
    return readout_dict

def readout(insertion_dict, chr_filter, min_matches=1):
    with open("magnify_readout.txt", "w") as file:
        file.write("Insertion Sites Found:\n")
        for read_chromosome, alignments in insertion_dict.items():
            file.write(read_chromosome+":\n")
            for align_chr, sites in alignments.items():
                if align_chr in chr_filter:
                    continue
                else:
                    file.write('\t'+align_chr+":\n")
                    for site, repetitions in sites.items():
                        if repetitions >= min_matches:
                            if str(site)[0] == '-':
                                file.write('\t'+'\t'+str(site)[1:]+" (Chimeric/Split Read): "+str(repetitions)+" matched\n")
                            else:
                                file.write('\t'+'\t'+str(site)+": "+str(repetitions)+" matched\n")
                        else:
                            continue

def regenerate_files(read_source,chr_list):
    source_name = f'magnify_{read_source}_markers_diff_chr.sam'
    source_reads = f'magnify_{read_source}_markers.sam'
    os.system(f"samtools view -H {source_reads} >> magnify_regenerated_temp.sam")
    with open(source_name) as sam:
        read_names = "|".join([SamAlignment(alignment).QNAME for alignment in sam])
        os.system(f"grep -E '{read_names}' {source_reads} >> magnify_regenerated_temp.sam")
    os.system('samtools sort magnify_regenerated_temp.sam > magnify_regenerated_temp.bam')
    os.system('samtools index magnify_regenerated_temp.bam')
    for chr in chr_list:
        os.system(f'samtools view -h magnify_regenerated_temp.bam {chr} > magnify_{chr}.sam')

#data = group('magnify_test_markers_diff_chr.sam')
#insertions = compress(data,granularity=0)
#readout(insertions,'',min_matches=1)
#regenerate_files('test',chr_list) #magnify_test_markers_diff_chr.sam magnify_test_markers.sam