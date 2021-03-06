#!/usr/bin/env python3
import sys
import os
cwd = os.getcwd() + "/automation"
sys.path.insert(0, cwd)
from automation_script import *
forwreads = []
revreads = []

#program folder locations
megahit_folder = '' #megahit folder
abyss_folder = '' #abyss folder
spades_folder = '' #spades folder
transrate_folder = '' #transrate folder
trimmomatic_folder = '' #trimmomatic folder
trinity_folder = '' #trinity location
salmon_folder = '' #salmon location
kallisto_folder = '' #kallisto location
bowtie_bin = '' #bowtie2 binary
blast_folder = '' #ncbi blast location
#---
custom_blast_folder = '' #custom blast library location
blast_name = '' #name of custom blast library; name in db folder
evalue = "1e-10" #blast e-value
trim_length = 300 #base pair length used in trim_contigs()

#change line below to pick which modules to run
def main():
    manual_input()
    trimmomatic()
    megahit()
    abyss()
    spades()
    transrate()
    trim_contigs()
    blastn()
    blastp()
    blastx()
    pull_matches() #use only if pull_matches_fast doesn't work
    pull_matches_fast()
    rsem()
    salmon()
    kallisto()

if __name__ == '__main__':
    main()


#input file location
input_srx = "path/to/file"
