'''
Functions to add include:
    - if foldername.txt file is not present skip makefolders() and fastqdump() and ask for the forward and reverse reads (so script can be run on individual data sets not published to NCBI)
    - run through all 3 assemblers in one go (Megahit, abyss, Spades)
    - configuration file for binary file locations
    - rework subprocess functions
    - run the pull seq python script within this scripts
    - seperate out each function of this file into their own files for better mobularity???
    - make sure it does not overwrite previously existing assemblies, just skips them
'''


import os
import shutil
import glob
import subprocess
from time import sleep
import sys
import pipeline_conf

# ########## User input
# input_srx = "/home/litoria/Documents/test/foldernames.txt"	#name of file with all SRX numbers and sample names
# spades_bin = "~/Assembly_Tools/SPAdes-3.12.0/bin/spades.py" #location of spades binary on system
# custom_location = "/home/litoria/Assembly_Tools/Uniprot_library" 	#location of custom blast library
# blast_name = "uniprot_db"		#name of custom blast library
# transrate_bin = "/home/litoria/Assembly_Tools/transrate-1.0.3-linux-x86_64/transrate"     #location of transrate binary on system
# megahit_bin = "megahit"
#
#
# blast_bin = "/home/litoria/NCBI_Tools/ncbi-blast-2.7.1+/bin/blastp"
# working_dir = "/home/litoria/Documents/test"
# # assembly_dir = f"{srx}.megahit_asm"
# # output_file = f"{subdir}_megahit.table"
# contig_file = []

#change names of all imported variables
if pipeline_conf.forwreads:
    forwreads = pipeline_conf.forwreads

if pipeline_conf.revreads:
    revreads = pipeline_conf.revreads

if pipeline_conf.megahit_bin:
    megahit_bin = pipeline_conf.megahit_bin

if pipeline_conf.abyss_bin:
    abyss_bin = pipeline_conf.abyss_bin

if pipeline_conf.spades_bin:
    spades_bin = pipeline_conf.abyss_bin

if pipeline_conf.transrate_bin:
    transrate_bin = pipeline_conf.transrate_bin

coreamount = int(os.cpu_count())
########## Functions


#make folders based on strings in input file
def makefolders(input_srx):
    global name
    with open(input_srx) as x:
        for line in x:
            name = line.strip()
            try:
                if not os.path.exists(os.getcwd() + "/" + name):
                    os.makedirs(os.getcwd() + "/" + name)
            except OSError:
                print("Error Making Directory" + name)

'''
Potentially add function here to grab SRX number by grabbing the last section of the directory being worked in
'''

#run fastqdump on specific srx number
def fastqdump(srx):
    global forwreads, revreads
    subprocess.call(["fastq-dump", "--dumpbase", "--defline-seq", "@$sn[_$rn]/$ri", "--split-files", srx])
    forwreads = f"{srx}_1.fastq"
    revreads = f"{srx}_2.fastq"

#manually select input files
def manual_input():
    global forwreads, revreads
    # forwreads = input("Drag forward reads here and press enter")
    # revreads = input("Drag reverse reads here and press enter")
    forwreads = sys.argv[1]
    revreads = sys.argv[2]


#run megahit in single or paired end mode depending on amount of input files found
def megahit1():     #use when manually inputing files
    global megahit_bin, forwreads, revreads
    if len(glob.glob("*.megahit_asm")) >= 1:
        print("Megahit assembly folder already present.\nSkipping Megahit Assembly")
        sleep(5)
    else:
        if not revreads:
            print("Running Megahit in single-end mode")
            sleep(5)
            subprocess.run(f"{megahit_bin} -r {forwreads} -o output.megahit_asm", shell=True)
        else:
            print("Running Megahit in paired-end mode")
            sleep(5)
            subprocess.run(f"{megahit_bin} -1 {forwreads} -2 {revreads} -o output.megahit_asm", shell=True)



def megahit():      #use when input files need to be automatically detected
    if len(glob.glob("*.megahit_asm")) >= 1:
        print("Megahit assembly folder already present")
        sleep(5)
    elif len(glob.glob("*.fastq")) == 2:
        subprocess.run(f"{megahit_bin} -1 {forwreads} -2 {revreads} -o {forwreads}.megahit_asm", shell=True)
    elif len(glob.glob("*.fastq")) == 1:
        subprocess.run(f"{megahit_bin} -r {forwreads} -o {forwreads}.megahit_asm", shell=True)
    else:
        print("Error: Input files were not found")

def abyss():
    global abyss_bin, forwreads, revreads
    if len(glob.glob("abyss")) >= 1:
        print("Abyss assembly folder already present.\nSkipping Abyss assembly")
        sleep(5)
    else:
        cur_dir = os.getcwd()
        os.mkdir(os.getcwd() + "/abyss")
        os.chdir(os.getcwd() + "/abyss")
        subprocess.run(f"{abyss_bin} j={coreamount} k=99 name=abyss_run in='{forwreads} {revreads}'", shell=True)
        os.chdir(cur_dir)

def spades():
    global spades_bin, forwreads, revreads
    subprocess.run(f"{spades_bin} -k 127 -1 {forwreads} -2 {revreads} -o spades_assembly", shell=True)


#combine all contig outputs with transrate
def transrate():
    global transrate_bin
    #find megahit contig file
    megahit_contig = (os.getcwd() + "/output.megahit_asm/final.contigs.fa")

    #find abyss contig file
    abyss_contig = (os.getcwd() + "/abyss/abyss_run-contigs.fa")

    #find spades contig file
    spades_contig = (os.getcwd() + "/spade_assembly/contigs.fasta")

    #run transrate to combine all contig files
    subprocess.run(f"{transrate_bin} --assembly {megahit_contig}, {abyss_contig}, {spades_contig} --merge-assemblies mergedassemblies", shell=True)


#copy custom blast library into working directory
def customblast():
    global custom_location, blast_name, contig_file
    try:
        contig_file = input("Drag contig file here and press enter:")
        custom_location = os.path.dirname(input("Drag custom blast library here and press enter:"))
        blast_name = input("Type name of custom blast library:")
        shutil.copy(contig_file, custom_location)

    # try:
    #     if not os.path.exists(os.getcwd() + "/blast_annotation"):
    #         os.mkdir(os.getcwd() + "/blast_annotation")
    #         custom_location = input("Drag custom blast library here and press enter")
    #         blast_name = [os.path.abspath(input("Type name of custom blast library:"))]
    #         shutil.copytree(custom_location, os.getcwd() + "/blast_annotation")
    except OSError:
        print("Blast directory could not be copied or already exists.")




#copy final.contigs.fa into customblast folder and run annotation
def annotation():
    global custom_location, blast_name, contig_file
    if not contig_file:
        contig_file = input("Drag contig file here and press enter")
    blast_code = [f"{blast_bin} -query {contig_file} -db {blast_name} -evalue 0.01 -max_target_seqs 1 -outfmt '7 std qseqid stitle sscinames staxids' -out output_blast.table -num_threads 12"]
    print(f"Assembly folder found: running annotation with {blast_name} library")

    current_dir = os.getcwd()
    os.chdir(custom_location)


    print(blast_code)
    subprocess.run(f"{blast_bin} -query {contig_file} -db {blast_name} -evalue 0.01 -max_target_seqs 1 -outfmt '7 std qseqid stitle sscinames staxids' -out output_blast.table -num_threads 12", shell=True)
    os.chdir(current_dir)


#create final output folder and move output files into this folder
def finalfolder():
    output_file = f"{subdir}_megahit.table"
    os.mkdir(subdir)
    os.chdir(subdir)
    shutil.copy(custom_location + "/" + filename, ".")
    shutil.copy(os.getcwd() + "/" + assembly_dir + filename, "." + subdir)


########## End of Functions


#make folders based on the lines within a textfile
# os.chdir(working_dir)
# makefolders(input_srx)
#
# subdirs = next(os.walk("."))[1]
# print(subdirs)												#remove after testing
# for subdir in subdirs:
# 	#grab srx from folder name
# 	print(subdir)
# 	srx = subdir.split("_")[2]
# 	print(srx)
# 	#open folder
# 	os.chdir(os.getcwd() + "/" + subdir)
# 	print(os.getcwd())										#remove after testing
# 	#run fastqdump
# 	fastqdump(srx)
# 	#run megahit
# 	megahit()
# 	#prepare folders for blast blast_annotation
# 	# customblast()
# 	#run annotation on custom blast library
# 	annotation()
# 	#make additional folder and copy final files into it
# 	finalfolder()
# 	#leave current folder to continue loop
# 	os.chdir("..")
# 	os.chdir("..")
