import os
import shutil
import glob
import subprocess
from time import sleep
import sys
import seqFrog_conf
import ntpath
import pandas as pd
import numpy as np
import fileinput
from pull_seqs_from_assembly import get_seqs
from multi import compare, get_seqs_fast
from filter_rsem import filter_rsem
from filter_length import filter_length



#change names of all imported variables
if seqFrog_conf.forwreads:
    forwreads = seqFrog_conf.forwreads

if seqFrog_conf.revreads:
    revreads = seqFrog_conf.revreads

if seqFrog_conf.trimmomatic_folder:
    trimmomatic_folder = seqFrog_conf.trimmomatic_folder

if seqFrog_conf.megahit_folder:
    megahit_folder = seqFrog_conf.megahit_folder

if seqFrog_conf.abyss_folder:
    abyss_folder = seqFrog_conf.abyss_folder

if seqFrog_conf.spades_folder:
    spades_folder = seqFrog_conf.spades_folder

if seqFrog_conf.transrate_folder:
    transrate_folder = seqFrog_conf.transrate_folder

if seqFrog_conf.custom_blast_folder:
    custom_blast_folder = seqFrog_conf.custom_blast_folder

if seqFrog_conf.blast_name:
    blast_name = seqFrog_conf.blast_name

if seqFrog_conf.blast_folder:
    blast_folder = seqFrog_conf.blast_folder

if seqFrog_conf.trinity_folder:
    trinity_folder = seqFrog_conf.trinity_folder

if seqFrog_conf.salmon_folder:
    salmon_folder = seqFrog_conf.salmon_folder

if seqFrog_conf.kallisto_folder:
    kallisto_folder = seqFrog_conf.kallisto_folder

if seqFrog_conf.bowtie_bin:
    bowtie_bin = seqFrog_conf.bowtie_bin

if seqFrog_conf.evalue:
    evalue = seqFrog_conf.evalue

if seqFrog_conf.trim_length:
    trim_length = seqFrog_conf.trim_length

#grab number of cpu cores for later use
coreamount = int(os.cpu_count())

contig_file = []
########## Functions

#manually select input files
def manual_input():
    global forwreads, revreads, contig_file, species_name, tissue_type, blast_file, one, two

    #grab arguments from console and pass them to python script
    import argparse
    parser = argparse.ArgumentParser()

    ####flags for the required inputs
    req_grp = parser.add_argument_group(title='required arguments')
    req_grp.add_argument("-f", "--forward-reads", dest = "forwreads", required=True, help="Forward Reads")
    req_grp.add_argument("-r", "--reverse-reads", dest = "revreads", required=True, help="Reverse Reads")
    ####end of flags for the required inputs

    #optional inputs
    parser.add_argument("-c", "--contigs", dest = "contig_file", help = "Contig file used if skipping assemblies", default = "")
    parser.add_argument("-sp", "--species", dest = "species_name", help = "Species name for organism being analysed (used for naming of final output files)", default = "")
    parser.add_argument("-t", "--tissue", dest = "tissue_type", help = "Tissue type for organism being analysed (used for naming of final output files)", default = "")
    parser.add_argument("-b", "-blast-file", dest = "blast_file", help = "Blast file for comparison", default = "")

    args = parser.parse_args()

    #assign variables from command line arguments
    forwreads = args.forwreads
    revreads = args.revreads
    contig_file = args.contig_file
    species_name = args.species_name
    tissue_type = args.tissue_type
    blast_file = args.blast_file

#Print which files were manually set
    print("\n")
    print(f"Forward reads:\n{forwreads}\n")
    print(f"Reverse reads:\n{revreads}\n")
    if contig_file:
        print(f"Contig file:\n{contig_file}\n")
    if blast_file:
        print(f"Blast file:\n{blast_file}\n")
    if species_name:
        one = f"[{species_name}]"
        print(f"Species name:\n{species_name}\n")
    else:
        one = ""
        print("Species not specified\n")

    if tissue_type:
        two = f"[{tissue_type}]"
        print(f"Tissue type:\n{tissue_type}\n")
    else:
        two = ""
        print("Tissue type not specified")

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
    global forwreads, revreads, trimmomatic_folder
    subprocess.run(f"java -jar {trimmomatic_folder}/trimmomatic-0.39.jar PE -threads {coreamount} -phred33 {forwreads} {revreads} forward_paired.fastq.gz forward_unpaired.fastq.gz reverse_paired.fastq.gz reverse_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36", shell=True)
    forwreads = (os.getcwd() +"/forward_paired.fastq.gz")
    revreads = (os.getcwd() + "/reverse_paired.fastq.gz")

def fastqdump(srx):
    global forwreads, revreads
    subprocess.call(["fastq-dump", "--dumpbase", "--defline-seq", "@$sn[_$rn]/$ri", "--split-files", srx])
    forwreads = f"{srx}_1.fastq"
    revreads = f"{srx}_2.fastq"


#run megahit in single or paired end mode depending on amount of input files found
def megahit():     #use when manually inputing files
    global megahit_folder, forwreads, revreads, contig_file
    if len(glob.glob("megahit_assembly")) >= 1:
        print("Megahit assembly folder already present.\nSkipping Megahit Assembly")
        sleep(5)
    else:
        if not revreads:
            print("Running Megahit in single-end mode")
            sleep(5)
            subprocess.run(f"{megahit_folder}/bin/megahit -r {forwreads} -o megahit_assembly", shell=True)
        else:
            print("Running Megahit in paired-end mode")
            sleep(5)
            subprocess.run(f"{megahit_folder}/bin/megahit -1 {forwreads} -2 {revreads} -o megahit_assembly", shell=True)

        contig_file = (os.getcwd() + "/megahit_assembly/final.contigs.fa")


def megahit1():      #use when input files need to be automatically detected NOT working
    '''
    add section to automatically pick forward and reverse revreads
    -if 1.fasta or 1.fa or 1.fa.gz then forward reads, same for reverse
    '''
    if len(glob.glob("megahit_assembly")) >= 1:
        print("Megahit assembly folder already present")
        sleep(5)
    elif len(glob.glob("*.fastq")) == 2:
        subprocess.run(f"{megahit_folder} -1 {forwreads} -2 {revreads} -o megahit_assembly", shell=True)
    elif len(glob.glob("*.fastq")) == 1:
        subprocess.run(f"{megahit_folder} -r {forwreads} -o megahit_assembly", shell=True)
    else:
        print("Error: Input files were not found")


def abyss():
    global abyss_folder, forwreads, revreads
    if len(glob.glob("abyss_assembly")) >= 1:
        print("Abyss assembly folder already present.\nSkipping Abyss assembly")
        sleep(5)
    else:
        cur_dir = os.getcwd()
        os.mkdir(os.getcwd() + "/abyss_assembly")
        os.chdir(os.getcwd() + "/abyss_assembly")
        subprocess.run(f"{abyss_folder}/bin/abyss-pe j={coreamount} k=59 name=abyss_run in='{forwreads} {revreads}'", shell=True)
        os.chdir(cur_dir)
    contig_file = (os.getcwd() + "/abyss_assembly/abyss_run-contigs.fa")



def spades():
    global spades_folder, forwreads, revreads
    if len(glob.glob("spades_assembly")) >= 1:
        print("Spades assembly folder already present.\nSkipping Spades assembly")
        sleep(5)

    else:
        subprocess.run(f"{spades_folder}/bin/spades.py --only-assembler -k 59 -1 {forwreads} -2 {revreads} -o spades_assembly", shell=True)
    contig_file = (os.getcwd() + "/spades_assembly/contigs.fasta")



#combine all contig outputs with transrate
def transrate():
    global transrate_folder, contig_file

    shutil.copy(os.getcwd() + "/megahit_assembly/final.contigs.fa", os.getcwd() + "/megahit_assembly/[backup]final.contigs.fa")
    with open(os.getcwd() + "/megahit_assembly/final.contigs.fa", 'r') as f:
        data = f.read()
    data = data.replace(" ","_")
    data = data.replace(",","_")
    with open(os.getcwd() + "/megahit_assembly/final.contigs.fa", 'w') as g:
        g.write(data)

    shutil.copy(os.getcwd() + "/abyss_assembly/abyss_run-contigs.fa", os.getcwd() + "/abyss_assembly/[backup]abyss_run-contigs.fa")
    with open(os.getcwd() + "/abyss_assembly/abyss_run-contigs.fa", 'r') as f:
        data = f.read()
    data = data.replace(" ","_")
    data = data.replace(",","_")
    with open(os.getcwd() + "/abyss_assembly/abyss_run-contigs.fa", 'w') as g:
        g.write(data)

    shutil.copy(os.getcwd() + "/spades_assembly/contigs.fasta", os.getcwd() + "/spades_assembly/[backup]contigs.fasta")
    with open(os.getcwd() + "/spades_assembly/contigs.fasta", 'r') as f:
        data = f.read()
    data = data.replace(" ","_")
    data = data.replace(",","_")
    with open(os.getcwd() + "/spades_assembly/contigs.fasta", 'w') as g:
        g.write(data)

    #find megahit contig file
    megahit_contig = (os.getcwd() + "/megahit_assembly/final.contigs.fa")

    #find abyss contig file
    abyss_contig = (os.getcwd() + "/abyss_assembly/abyss_run-contigs.fa")

    #find spades contig file
    spades_contig = (os.getcwd() + "/spades_assembly/contigs.fasta")

    if len(glob.glob("transrate_results")) >= 1:
        print("Transrate folder already present.\nSkipping Transrate")
        contig_file = (os.getcwd() + "/transrate_results/merged_assemblies")
        sleep(5)
    else:
        #run transrate to combine all contig files
        subprocess.run(f"{transrate_folder}/transrate --assembly {megahit_contig},{abyss_contig},{spades_contig} --merge-assemblies merged_assemblies", shell=True)

        contig_file = (os.getcwd() + "/transrate_results/merged_assemblies")

def trim_contigs():
    global contig_file
    filter_length(contig_file, trim_length)
    os.rename("good.fasta", "good_merged_assemblies")
    contig_file = (os.getcwd() + "/good_merged_assemblies")


#run annotation with database set in config file
def blastn():
    global contig_file, custom_blast_folder, blast_name, species_name, tissue_type, blast_file, one, two
    '''
    add section to do annotation from individual assemblies if merge was not needed (smaller data sets)
    '''

    if not contig_file:
        print("No contig file was found.")
        return #stops program if contig_file is empty/not set
    else:
        print(f"Assembly folder found: running annotation on {contig_file} with {blast_name} library")

        wd = os.getcwd()
        os.chdir(custom_blast_folder)
        filename = ntpath.basename(f"{contig_file}")
        print(filename)

        subprocess.run(f"{blast_folder}/bin/blastn -query {contig_file} -db {blast_name} -evalue {evalue} -max_target_seqs 1 -outfmt '7 std qseqid stitle sscinames staxids' -out {one}{two}[{blast_name}]{filename}_blast.table -num_threads 12", shell=True)
        blast_file = (f"{one}{two}[{blast_name}]{filename}_blast.table")
        shutil.move(blast_file, wd)
        os.chdir(wd)

def blastp():
    global contig_file, custom_blast_folder, blast_name, species_name, tissue_type, blast_file, one, two
    '''
    add section to do annotation from individual assemblies if merge was not needed (smaller data sets)
    '''

    if not contig_file:
        print("No contig file was found.")
        return #stops program if contig_file is empty/not set
    else:
        print(f"Assembly folder found: running annotation on {contig_file} with {blast_name} library")

        wd = os.getcwd()
        os.chdir(custom_blast_folder)
        filename = ntpath.basename(f"{contig_file}")
        print(filename)

        subprocess.run(f"{blast_folder}/bin/blastp -query {contig_file} -db {blast_name} -evalue {evalue} -max_target_seqs 1 -outfmt '7 std qseqid stitle sscinames staxids' -out {one}{two}[{blast_name}]{filename}_blast.table -num_threads 12", shell=True)
        blast_file = (f"{one}{two}[{blast_name}]{filename}_blast.table")
        shutil.move(blast_file, wd)
        os.chdir(wd)


def blastx():
    global contig_file, custom_blast_folder, blast_name, species_name, tissue_type, blast_file, one, two

    if not contig_file:
        print("No contig file was found.")
        return #stops program if contig_file is empty/not set
    else:
        print(f"Assembly folder found: running BLASTx annotation on {contig_file} with {blast_name} library")

        wd = os.getcwd()
        os.chdir(custom_blast_folder)
        filename = ntpath.basename(f"{contig_file}")
        print(filename)

        subprocess.run(f"{blast_folder}/bin/blastx -query {contig_file} -db {blast_name} -evalue {evalue} -max_target_seqs 1 -outfmt '7 std qseqid stitle sscinames staxids' -out {one}{two}[{blast_name}]{filename}_blast.table -num_threads {coreamount}", shell=True)
        blast_file = (f"{one}{two}[{blast_name}]{filename}_blast.table")
        shutil.move(blast_file, wd)
        os.chdir(wd)


def pull_matches():
    global contig_file, blast_file
    full_name = f"{blast_file}_matches.fasta"

    comparecolumn = "1"
    keepcolumn = "1,2,3,4,5"
    get_seqs(contig_file, blast_file, full_name, [int(x) for x in comparecolumn.split(",")], [int(x) for x in keepcolumn.split(",")])


def pull_matches_fast():
    global contig_file, blast_file
    full_name = f"{blast_file}_matches.fasta"

    comparecolumn = "1"
    keepcolumn = "1,2,3,4,5"
    get_seqs_fast(contig_file, blast_file, full_name, [int(x) for x in comparecolumn], [int(x) for x in keepcolumn.split(",")])

def rsem():
    with fileinput.FileInput(f"{blast_file}_matches.fasta", inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(" ", "_").replace("\t", "__"), end='')
    os.rename(f"{blast_file}_matches.fasta", "matches.fasta")

    subprocess.run(f"{trinity_folder}/util/align_and_estimate_abundance.pl --transcripts matches.fasta --seqType fq --left {forwreads} --right {revreads} --est_method RSEM --aln_method {bowtie_bin} --prep_reference --output_dir rsem_results", shell=True)
    rsem_file = (os.getcwd() + "/rsem_results/RSEM.genes.results")
    filter_rsem(rsem_file)
    os.rename("matches.fasta", f"{blast_file}_matches.fasta")


def salmon():
    with fileinput.FileInput(f"{blast_file}_matches.fasta", inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(" ", "_").replace("\t", "__"), end='')
    os.rename(f"{blast_file}_matches.fasta", "matches.fasta")

    index = f"{species_name}{tissue_type}_index"

    subprocess.run(f"{salmon_folder}/bin/salmon index -t matches.fasta -i {index}", shell=True)
    subprocess.run(f"{salmon_folder}/bin/salmon quant -i {index} -l A -1 {forwreads} -2 {revreads} -p {coreamount} -o quants/{species_name}{tissue_type}", shell=True)
    os.rename("matches.fasta", f"{blast_file}_matches.fasta")


def kallisto():
    with fileinput.FileInput(f"{blast_file}_matches.fasta", inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(" ", "_").replace("\t", "__"), end='')
    os.rename(f"{blast_file}_matches.fasta", "matches.fasta")

    cwd = os.getcwd()
    os.mkdir(os.getcwd() + "/kallisto")
    os.chdir(os.getcwd() + "/kallisto")

    subprocess.run(f"{kallisto_folder}/kallisto index -i transcripts.idx '{cwd}/matches.fasta'", shell=True)
    subprocess.run(f"{kallisto_folder}/kallisto quant -i transcripts.idx -o output -b 100 {forwreads} {revreads}", shell=True)

    os.chdir(cwd)
    os.rename("matches.fasta", f"{blast_file}_matches.fasta")


########## End of Functions
