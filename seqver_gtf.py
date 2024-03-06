from seqver_genomeupdate import *

# --------------------------------------------------------------------
# | Guide to cases (a-f) as described in this program                |
# | ~ denotes gtf feature, + denotes deletion, - denotes other bases |
# --------------------------------------------------------------------
#
# (a): edit strictly before feature - this needs to be shifted depending on the insertion
# +++++
# ------~~~~~~~~~~------
#
# (b): edit overlapping with feature on the left - this needs to trim the left side of the insertion and then shift the feature depending on the insertion
#     ++++++
# ------~~~~~~~~~~------
#
# (c): edit enveloping entire feature - this needs to delete the feature
#      ++++++++++++
# ------~~~~~~~~~~------
#
# (d): edit completely internal to feature - this needs to shift the end of the feature depending on the length of the insertion
#         +++++
# ------~~~~~~~~~~------
#
# (e): edit overlapping with feature on the right - this needs to trim the right side of the insertion but not shift the feature
#               +++++
# ------~~~~~~~~~~------
#
# (f): edit strictly after feature - this needs to leave the feature untouched since everything happens after it
#                  +++++
# ------~~~~~~~~~~------

def gtfEdit(gtf_file, commands, temp_folder, folder, unsorted_output, sorted_output):
    with open(gtf_file,"r") as gtf:
        with open(f"{temp_folder}/{unsorted_output}","w") as out:
            new_gtfs = set()
            for line in gtf:
                if not line.startswith("#"):
                    deleteFlag = False
                    fields = line.split("\t")
                    gtf_chr, gtf_start, gtf_end = fields[0], int(fields[3]), int(fields[4])
                    for command in commands:
                        command = str(command)
                        command_fields, insertion_len, insertion_name = command.split("\t")[0], len(command.split("\t")[1]), str(command.split("\t")[2])
                        cmd_chr, cmd_start, cmd_end = command_fields.split(":")[0], int(command_fields.split(":")[1].split("-")[0]), int(command_fields.split(":")[1].split("-")[1])
                        shift = insertion_len - (int(cmd_end) - int(cmd_start))
                        if gtf_chr == cmd_chr:
                            if gtf_start > cmd_start: #have to be unequal otherwise you could get internal insertions with (c) which we don't want
                            #the -1 is because cmd coordinates are 0-based, gtf coordinates are 1-based
                                if gtf_start <= cmd_end: #equal because command endings are included so if they're the same it means we want to trim tart by 1
                                    if gtf_end <= cmd_end: #equal since endpoints are included, so if the two ends match then it's still a complete deletion
                                        #(c) completely deletes a feature
                                        deleteFlag = True 
                                    else: 
                                        #(b) trim start of the feature
                                        trimming_factor = gtf_start - (cmd_start+1) #+1 here is due to the fact that the starting coordinate isn't included in commands
                                        gtf_start -= trimming_factor
                                        gtf_start += insertion_len
                                        gtf_end += shift
                                elif gtf_start > cmd_end: #unequal because equality here means a deletion of the first coordinate so it's not only shifting.
                                    #(a) just shift depending on the insertion,
                                    gtf_start += shift
                                    gtf_end += shift
                            elif gtf_start <= cmd_start: #equal because starting coords aren't included in commands so if the two are equal, the insertion starts at start+1 so it's internal
                                if gtf_end < cmd_end: #unequal because endpoints are included so equality means we want (d), an edit internal to the feature.
                                    if gtf_end <= cmd_start: #equal since commands don't include start point, so if they're equal the insertion starts after the feature
                                        #(f) don't do anything since all changes happen after the gtf
                                        pass
                                    elif gtf_end > cmd_start: #unequal since commands don't include start point, so if they're equal we want (f), all of the command is external
                                        #(e) trim end of feature 
                                        trimming_factor = gtf_end - cmd_start
                                        gtf_end -= trimming_factor
                                        #NOTE: no shift, since any insertions would happen after the end of the feature and therefore the feature does not need to shift to account for them
                                elif gtf_end >= cmd_end: 
                                    #(d) edits internal to the feature
                                    trimming_factor = cmd_end-cmd_start
                                    gtf_end -= trimming_factor
                                    gtf_end += insertion_len
                    fields[3], fields[4] = str(gtf_start), str(gtf_end) #no empty lines in gtf file are allowed
                    fields[-1] = str(fields[-1].strip())+"\n" #makes sure the gtf line ends in a newline to ensure new features start on a new line as intended
                    if not deleteFlag:
                        out.write("\t".join(fields))
                    if insertion_len > 0:
                        new_feature_fields = fields
                        new_feature_fields[1], new_feature_fields[2], new_feature_fields[3], new_feature_fields[4], new_feature_fields[8] = "seqverify", "gene", str(cmd_start+1), str(cmd_end+shift), insertion_name
                        new_feature_fields[-1] = str(new_feature_fields[-1].strip())+"\n"
                        new_line = "\t".join(new_feature_fields)
                        feature_unique = "\t".join(new_feature_fields[0:4])
                        if feature_unique not in new_gtfs:
                            out.write(new_line)
                            new_gtfs.add(feature_unique)
        os.system(f'''awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k4,4n -k5,5n"}}' {temp_folder}/{unsorted_output} > {folder}/{sorted_output}''') #~20min to sort on my laptop at 15GB ram, 2.3GHz