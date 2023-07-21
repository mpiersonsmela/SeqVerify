import math as m
import os

class editArgument:
    def __init__(self, argument):
        fields = argument.split("\t")
        site = fields[0].split(":")
        self.chr = site[0]
        self.coords = [int(i) for i in site[1].split("-")] #translates between SAM coords (1-based) and Python coords (0-based), assumes input in SAM coords
        self.sequence = fields[1].strip()
    
    def __str__(self):
        return f"{self.chr}:{self.coords[0]}-{self.coords[1]}\t{self.sequence}"
    
    def iterable(self):
        return [self.chr,self.coords[0],self.coords[1],self.sequence]

def selectFromPosition(fasta, chr): #takes the relevant sequence info from an arbitrary FASTA file
    with open(f"{fasta}","r") as reference:
        relevant, startSeen = [], False
        for line in reference:
            if startSeen:
                if not line.startswith(">"):
                    relevant.append(line)
                else:
                    break
            else:
                if line.startswith(f">{chr} ") or line.startswith(f">{chr}\n"):
                    startSeen = True
                else:
                    continue
    max_len = max([len(i) for i in relevant])-1
    relevant = [i.replace("\n","") for i in relevant]
    return [relevant, max_len]

def genomeSplitter(genome,chr,folder):
    with open(f"{folder}/{chr}.fa","w+") as chr_seq:
        chr_seq.write(f">{chr}\n")
        chr_seq.writelines([i+"\n" for i in selectFromPosition(genome,chr)[0]])

def chrCommand(command,folder):
    sequence = f"{command.chr}.fa"
    position = selectFromPosition(os.getcwd()+f"/{folder}/"+sequence,command.chr) 
    replacement = "".join(position[0])
    replacement = replacement[0:command.coords[0]]+command.sequence+replacement[command.coords[1]:]
    replacement = [replacement[(i)*position[1]:(i+1)*position[1]] for i in range(m.ceil(len(replacement)/position[1]))]
    with open(f"{folder}/{sequence}","w+") as reference:
        reference.write(f">{command.chr}\n")
        reference.writelines([i+"\n" for i in replacement])

def commandHandler(genome, commandfile, folder, newName): #assumes that within each chromosome, commands are non-overlapping and unique 
    with open(f"{commandfile}","r") as original:
        commands = original.readlines()
        commandfile = commandfile.split("/")[-1]
        arguments = [editArgument(command) for command in commands]
        chrsAffected = list(set([argument.chr for argument in arguments]))
        commandGroup = {chr: [argument for argument in arguments if argument.chr == chr] for chr in chrsAffected}
        commandGroup = {chr: sorted(commandGroup[chr], key= lambda x: x.iterable()[1]) for chr in commandGroup.keys()}
        for chr in commandGroup.keys():
            for i,cmd in enumerate(commandGroup[chr]):
                shift = (cmd.coords[0]-cmd.coords[1]) + len(cmd.sequence) 
                for cmdToEdit in commandGroup[chr][i+1:]:
                    cmdToEdit.coords = [coord+shift for coord in cmdToEdit.coords]
        totalCommands = [command for group in commandGroup.values() for command in group]

        for chr in chrsAffected:
            genomeSplitter(genome, chr, folder)

        for command in totalCommands:
            chrCommand(command, folder)

        with open(f'{folder}/{newName}',"w+") as newGenome:
            os.system(f"grep -F '>' {genome} > {folder}/{newName}")
            refChrs = [i.strip() for i in newGenome.readlines()]
            newGenome.seek(0)
            newGenome.truncate(0)
            for chr in refChrs:
                chrName = chr.split(" ")[0][1:]
                if chrName in chrsAffected:
                    with open(f'{folder}/{chrName}.fa',"r") as affected:
                        toWrite = affected.readlines()
                else:
                    name = f"{chr}"
                    toWrite = selectFromPosition(genome,chr[1:])[0]
                    toWrite.insert(0,name)
                    toWrite = [i+"\n" for i in toWrite]
                newGenome.writelines(toWrite)    
    return totalCommands