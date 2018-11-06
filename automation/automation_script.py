'''
Functions to add include:
    - if foldername.txt file is not present skip makefolders() and fastqdump() and ask for the forward and reverse reads (so script can be run on individual data sets not published to NCBI)
    - rework subprocess functions
    - run the pull seq python script within this scripts
    - seperate out each function of this file into their own files for better modularity???
'''


import os
import shutil
import glob
import subprocess
from time import sleep
import sys
import pipeline_conf
import test

#change names of all imported variables
# if pipeline_conf.forwreads:
#     forwreads = pipeline_conf.forwreads
#
# if pipeline_conf.revreads:
#     revreads = pipeline_conf.revreads
if test.forwreads:
    forwreads = test.forwreads

if test.revreads:
    revreads = test.revreads


if pipeline_conf.megahit_bin:
    megahit_bin = pipeline_conf.megahit_bin

if pipeline_conf.abyss_bin:
    abyss_bin = pipeline_conf.abyss_bin

if pipeline_conf.spades_bin:
    spades_bin = pipeline_conf.spades_bin

if pipeline_conf.transrate_bin:
    transrate_bin = pipeline_conf.transrate_bin

if pipeline_conf.custom_location:
    custom_location = pipeline_conf.custom_location

if pipeline_conf.blast_name:
    blast_name = pipeline_conf.blast_name

if pipeline_conf.blast_bin:
    blast_bin = pipeline_conf.blast_bin

if pipeline_conf.trimmomatic_jar:
    trimmomatic_jar = pipeline_conf.trimmomatic_jar

coreamount = int(os.cpu_count())
########## Functions

#manually select input files
def manual_input():
    global forwreads, revreads
    # forwreads = input("Drag forward reads here and press enter")
    # revreads = input("Drag reverse reads here and press enter")
    forwreads = sys.argv[1]
    revreads = sys.argv[2]
    print(f"Forward reads:\n{forwreads}\n")
    print(f"Reverse reads:\n{revreads}\n")
    sleep(5)


#make folders based on strings in input file; DOES NOT WORK
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
def trimmomatic():
    global forwreads, revreads, trimmomatic_jar
    subprocess.run(f"java -jar {trimmomatic_jar} PE -threads 10 -phred33 {forwreads} {revreads} forward_paired.fastq forward_unpaired.fastq reverse_paired.fastq reverse_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36", shell=True)


def fastqdump(srx):
    global forwreads, revreads
    subprocess.call(["fastq-dump", "--dumpbase", "--defline-seq", "@$sn[_$rn]/$ri", "--split-files", srx])
    forwreads = f"{srx}_1.fastq"
    revreads = f"{srx}_2.fastq"


#run megahit in single or paired end mode depending on amount of input files found
def megahit1():     #use when manually inputing files
    global megahit_bin, forwreads, revreads
    if len(glob.glob("megahit_assembly")) >= 1:
        print("Megahit assembly folder already present.\nSkipping Megahit Assembly")
        sleep(5)
    else:
        if not revreads:
            print("Running Megahit in single-end mode")
            sleep(5)
            subprocess.run(f"{megahit_bin} -r {forwreads} -o megahit_assembly", shell=True)
        else:
            print("Running Megahit in paired-end mode")
            sleep(5)
            subprocess.run(f"{megahit_bin} -1 {forwreads} -2 {revreads} -o megahit_assembly", shell=True)

        contig_file = (os.getcwd() + "/megahit_assembly/final.contigs.fa")
        return contig_file


def megahit():      #use when input files need to be automatically detected NOT working
    '''
    add section to automatically pick forward and reverse revreads
    -if 1.fasta or 1.fa or 1.fa.gz then forward reads, same for reverse
    '''
    if len(glob.glob("megahit_assembly")) >= 1:
        print("Megahit assembly folder already present")
        sleep(5)
    elif len(glob.glob("*.fastq")) == 2:
        subprocess.run(f"{megahit_bin} -1 {forwreads} -2 {revreads} -o megahit_assembly", shell=True)
    elif len(glob.glob("*.fastq")) == 1:
        subprocess.run(f"{megahit_bin} -r {forwreads} -o megahit_assembly", shell=True)
    else:
        print("Error: Input files were not found")


def abyss():
    global abyss_bin, forwreads, revreads
    if len(glob.glob("abyss_assembly")) >= 1:
        print("Abyss assembly folder already present.\nSkipping Abyss assembly")
        sleep(5)
    else:
        cur_dir = os.getcwd()
        os.mkdir(os.getcwd() + "/abyss_assembly")
        os.chdir(os.getcwd() + "/abyss_assembly")
        subprocess.run(f"{abyss_bin} j={coreamount} k=99 name=abyss_run in='{forwreads} {revreads}'", shell=True)
        os.chdir(cur_dir)


def spades():
    global spades_bin, forwreads, revreads
    if len(glob.glob("spades_assembly")) >= 1:
        print("Spades assembly folder already present.\nSkipping Spades assembly")
        sleep(5)

    else:
        subprocess.run(f"{spades_bin} --only-assembler -k 99 -1 {forwreads} -2 {revreads} -o spades_assembly", shell=True)


#combine all contig outputs with transrate
def transrate():
    global transrate_bin
    #find megahit contig file
    megahit_contig = (os.getcwd() + "/megahit_assembly/final.contigs.fa")

    #find abyss contig file
    abyss_contig = (os.getcwd() + "/abyss_assembly/abyss_run-contigs.fa")

    #find spades contig file
    spades_contig = (os.getcwd() + "/spade_assembly/contigs.fasta")

    #run transrate to combine all contig files
    subprocess.run(f"{transrate_bin} --assembly {megahit_contig},{abyss_contig},{spades_contig} --merge-assemblies merged_assemblies", shell=True)

    contig_file = (os.getcwd() + "/transrate_results/merged_assemblies")
    return contig_file #output of this function will be this variable


#run annotation with database set in config file
def annotation():
    '''
    add section to do annotation from individual assemblies if merge was not needed (smaller data sets)
    '''
    if not contig_file:
        contig_file = (os.getcwd() + "/megahit_assembly/final.contigs.fa")
    global custom_location, blast_name
    if not contig_file:
        print("No contig file was found.")
        return #stops program if contig_file is empty/not set
    print(f"Assembly folder found: running annotation with {blast_name} library")

    wd = os.getcwd()
    os.chdir(custom_location)

    subprocess.run(f"{blast_bin} -query {contig_file} -db {blast_name} -evalue 0.01 -max_target_seqs 1 -outfmt '7 std qseqid stitle sscinames staxids' -out {contig_file}_blast.table -num_threads 12", shell=True)
    shutil.copy(f"{contig_file}_blast.table", wd)
    os.chdir(wd)


'''
Add section to pick out matches from annotation
remove remove_spaces
run rsem
output expression files
'''

########## End of Functions
