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

import os
import subprocess

def gtfEdit(gtf_file, commands, temp_folder, folder, unsorted_output, sorted_output):
    print(f"Starting GTF editing process")
    print(f"Input GTF file: {gtf_file}")
    print(f"Number of commands: {len(commands)}")
    print(f"Temporary folder: {temp_folder}")
    print(f"Output folder: {folder}")
    
    if not os.path.exists(gtf_file):
        print(f"Error: Input GTF file does not exist: {gtf_file}")
        return False
    
    if not os.access(temp_folder, os.W_OK):
        print(f"Error: No write permission for temporary folder: {temp_folder}")
        return False
    
    if not os.access(folder, os.W_OK):
        print(f"Error: No write permission for output folder: {folder}")
        return False
    
    try:
        with open(gtf_file, "r") as gtf:
            with open(f"{temp_folder}/{unsorted_output}", "w") as out:
                new_gtfs = set()
                lines_processed = 0
                for line in gtf:
                    lines_processed += 1
                    if lines_processed % 100000 == 0:
                        print(f"Processed {lines_processed} lines")
                    
                    if not line.startswith("#"):
                        deleteFlag = False
                        fields = line.split("\t")
                        gtf_chr, gtf_start, gtf_end = fields[0], int(fields[3]), int(fields[4])
                        for command in commands:
                            command = str(command)
                            command_parts = command.split("\t")
                            if len(command_parts) != 3:
                                print(f"Warning: Skipping invalid command: {command}")
                                continue
                            command_fields, insertion_len, insertion_name = command_parts[0], len(command_parts[1]), str(command_parts[2])
                            cmd_parts = command_fields.split(":")
                            if len(cmd_parts) != 2 or "-" not in cmd_parts[1]:
                                print(f"Warning: Skipping command with invalid format: {command}")
                                continue
                            cmd_chr, cmd_coords = cmd_parts[0], cmd_parts[1].split("-")
                            if len(cmd_coords) != 2:
                                print(f"Warning: Skipping command with invalid coordinates: {command}")
                                continue
                            cmd_start, cmd_end = int(cmd_coords[0]), int(cmd_coords[1])
                            shift = insertion_len - (int(cmd_end) - int(cmd_start))
                            
                            # ... (rest of the logic remains the same)
                    
                    fields[3], fields[4] = str(gtf_start), str(gtf_end)
                    fields[-1] = str(fields[-1].strip()) + "\n"
                    if not deleteFlag:
                        out.write("\t".join(fields))
                    if insertion_len > 0:
                        new_feature_fields = fields.copy()
                        new_feature_fields[1], new_feature_fields[2], new_feature_fields[3], new_feature_fields[4], new_feature_fields[8] = "seqverify", "gene", str(cmd_start+1), str(cmd_end+shift), insertion_name
                        new_feature_fields[-1] = str(new_feature_fields[-1].strip()) + "\n"
                        new_line = "\t".join(new_feature_fields)
                        feature_unique = "\t".join(new_feature_fields[0:4])
                        if feature_unique not in new_gtfs:
                            out.write(new_line)
                            new_gtfs.add(feature_unique)
                
                print(f"Finished processing GTF. Total lines processed: {lines_processed}")
        
        print(f"Sorting GTF file")
        sort_command = f'''awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k4,4n -k5,5n"}}' {temp_folder}/{unsorted_output} > {folder}/{sorted_output}'''
        subprocess.run(sort_command, shell=True, check=True)
        print(f"Sorting complete. Output file: {folder}/{sorted_output}")
        
        if os.path.exists(f"{folder}/{sorted_output}"):
            print(f"GTF editing process completed successfully")
            return True
        else:
            print(f"Error: Sorted output file not created: {folder}/{sorted_output}")
            return False
    
    except Exception as e:
        print(f"Error during GTF editing: {str(e)}")
        return False