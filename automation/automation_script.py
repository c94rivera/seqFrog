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
import ntpath
sys.path.insert(0, '/home/litoria/Documents/git/transcriptome_assembly/post_processing/')
import pull_seqs_from_assembly



#change names of all imported variables
# if pipeline_conf.forwreads:
#     forwreads = pipeline_conf.forwreads
#
# if pipeline_conf.revreads:
#     revreads = pipeline_conf.revreads
if pipeline_conf.forwreads:
    forwreads = pipeline_conf.forwreads

if pipeline_conf.revreads:
    revreads = pipeline_conf.revreads


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
contig_file = []
########## Functions

#manually select input files
def manual_input():
    global forwreads, revreads, contig_file, species_name, tissue_type

    #grab arguments from console and pass them to python script
    import argparse
    parser = argparse.ArgumentParser()

    #argument tags
    # parser.add_argument("contigs", help = "Contigs Input file")
    # parser.add_argument("blast_hits", help = "Blast Input")

    ####flags for the required inputs
    req_grp = parser.add_argument_group(title='required arguments')
    req_grp.add_argument("-f", "--forward-reads", dest = "forwreads", required=True, help="Forward Reads")
    req_grp.add_argument("-r", "--reverse-reads", dest = "revreads", required=True, help="Reverse Reads")
    ####end of flags for the required inputs

    parser.add_argument("-c", "--contigs", dest = "contig_file", help = "Contig file used if skipping assemblies", default = "")
    parser.add_argument("-sp", "--species", dest = "species_name", help = "Species name for organism being analysed (used for naming of final output files)", default = "")
    parser.add_argument("-t", "--tissue", dest = "tissue_type", help = "Tissue type for organism being analysed (used for naming of final output files)", default = "")

    args = parser.parse_args()

    #assign variables from command line arguments
    forwreads = args.forwreads
    revreads = args.revreads
    contig_file = args.contig_file
    species_name = args.species_name
    tissue_type = args.tissue_type
    blast_file = []

    print("\n")
    print(f"Forward reads:\n{forwreads}\n")
    print(f"Reverse reads:\n{revreads}\n")
    print(f"Contig file:\n{contig_file}\n")
    print(f"Species name:\n{species_name}\n")
    print(f"Tissue type:\n{tissue_type}\n")

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
    subprocess.run(f"java -jar {trimmomatic_jar} PE -threads 10 -phred33 {forwreads} {revreads} forward_paired.fastq.gz forward_unpaired.fastq.gz reverse_paired.fastq.gz reverse_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36", shell=True)
    forwreads = (os.getcwd() +"/forward_paired.fastq.gz")
    revreads = (os.getcwd() + "/reverse_paired.fastq.gz")

def fastqdump(srx):
    global forwreads, revreads
    subprocess.call(["fastq-dump", "--dumpbase", "--defline-seq", "@$sn[_$rn]/$ri", "--split-files", srx])
    forwreads = f"{srx}_1.fastq"
    revreads = f"{srx}_2.fastq"


#run megahit in single or paired end mode depending on amount of input files found
def megahit1():     #use when manually inputing files
    global megahit_bin, forwreads, revreads, contig_file
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
        subprocess.run(f"{abyss_bin} j={coreamount} k=59 name=abyss_run in='{forwreads} {revreads}'", shell=True)
        os.chdir(cur_dir)


def spades():
    global spades_bin, forwreads, revreads
    if len(glob.glob("spades_assembly")) >= 1:
        print("Spades assembly folder already present.\nSkipping Spades assembly")
        sleep(5)

    else:
        subprocess.run(f"{spades_bin} --only-assembler -k 59 -1 {forwreads} -2 {revreads} -o spades_assembly", shell=True)


#combine all contig outputs with transrate
def transrate():
    global transrate_bin, contig_file
    #find megahit contig file
    megahit_contig = (os.getcwd() + "/megahit_assembly/final.contigs.fa")

    #find abyss contig file
    abyss_contig = (os.getcwd() + "/abyss_assembly/abyss_run-contigs.fa")

    #find spades contig file
    spades_contig = (os.getcwd() + "/spades_assembly/contigs.fasta")

    if len(glob.glob("transrate_results")) >= 1:
        print("Transrate folder already present.\nSkipping Transrate")
        sleep(5)
    else:
        #run transrate to combine all contig files
        subprocess.run(f"{transrate_bin} --assembly {megahit_contig},{abyss_contig},{spades_contig} --merge-assemblies merged_assemblies", shell=True)

        contig_file = (os.getcwd() + "/transrate_results/merged_assemblies")


#run annotation with database set in config file
def annotation():
    global contig_file, custom_location, blast_name, species_name, tissue_type, blast_file
    '''
    add section to do annotation from individual assemblies if merge was not needed (smaller data sets)
    '''

    if not contig_file: #why????
        contig_file = (os.getcwd() + "/megahit_assembly/final.contigs.fa")

    if not contig_file:
        print("No contig file was found.")
        return #stops program if contig_file is empty/not set
    else:
        print(f"Assembly folder found: running annotation on {contig_file} with {blast_name} library")

        wd = os.getcwd()
        os.chdir(custom_location)
        filename = ntpath.basename(f"{contig_file}")
        print(filename)

        if species_name:
            one = f"[{species_name}]"
            print(species_name)
        else:
            one = ""
            print("Species not specified")

        if tissue_type:
            two = f"[{tissue_type}]"
            print(tissue_type)
        else:
            two = ""
            print("Tissue type not specified")


        subprocess.run(f"{blast_bin} -query {contig_file} -db {blast_name} -evalue 0.01 -max_target_seqs 1 -outfmt '7 std qseqid stitle sscinames staxids' -out {one}{two}[{blast_name}]{filename}_blast.table -num_threads 12", shell=True)
        blast_file = (f"{one}{two}[{blast_name}]{filename}_blast.table")
        shutil.move(blast_file, wd)
        os.chdir(wd)


def pull_matches():
    full_name = (f"{one}{two}[{blast_name}]{filename}_blast.table")
    get_seqs(contig_file, blast_file, full_name, "1", "1,2,3,10,13")
    '''RUN PYTHON SCRIPT ADD LENGTH OF CONTIGS FOR ANALYSIS
    NAME MATCHES'''
#
#
#
#
# def rsem():
#
#
#
'''
Add section to pick out matches from annotation
remove remove_spaces
run rsem
output expression files
'''

########## End of Functions
