def mutation_logger(folder,output_file,VCF_file,min_intensity="MODERATE"):
    with open(f"{folder}/{VCF_file}","r") as vcf:
        intensities = ["MODIFIER","LOW","MODERATE","HIGH"]
        min_intensity_index = intensities.index(min_intensity)
        with open(f"{folder}/{output_file}","w+") as output:
            output.write("CHR:COORD\tVARIANT_TYPE\tEFFECT\tGENE\tREFSEQ_ID\tDNA_change\tPROTEIN_change\tHOMOZYGOUS\tLOF\n")
            for line in vcf:
                if not line.startswith("#"):
                    fields = line.split("\t")
                    coords,annotations,genotype = ":".join(fields[0:2]),fields[7].split(";"),fields[-1].split(":")[0]
                    homozygous = (len(set(genotype.split("/"))) == 1)
                    
                    relevant = [i for i in annotations if i.startswith("ANN") or i.startswith("LOF")]
                    annotations = relevant[0].split(",")
                    try:
                        LOF_gene = relevant[1].split("|")[1]
                        LOF_gene = True
                    except:
                        LOF_gene = False

                    for ann in annotations:
                        args = ann.split("|")

                        if args[2] not in intensities[0:min_intensity_index]:
                            id_str = coords,args[1],args[2],args[3],args[6],args[9], args[10], homozygous,LOF_gene
                            id_str = [str(i) for i in id_str]
                            output.write("\t".join(id_str)+"\n")