#!/usr/bin/env python3
import sys
import os
cwd = os.getcwd() + "/automation"
sys.path.insert(0, cwd)
from automation_script import *
forwreads = []
revreads = []

#program folder locations
megahit_folder = "/megahit" #megahit binary
abyss_folder = '/home/phyllobates/Assembly_Tools/abyss-2.1.5/ABYSS/abyss.cc' #abyss binary
spades_folder = '/home/phyllobates/Assembly_Tools/SPAdes-3.13.0-Linux/bin/spades.py' #spades binary
transrate_folder = '/home/phyllobates/Assembly_Tools/transrate-1.0.3-linux-x86_64/transrate' #transrate binary
trimmomatic_folder = '/home/phyllobates/Assembly_Tools/Trimmomatic-0.38/trimmomatic-0.38.jar' #trimmomatic binary
#rsem_loc should be the location of the "align_and_estimate_abundance.pl" script within your Trinity install
trinity_folder = '/home/litoria/Assembly_Tools/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl' #rsem location
salmon_folder = 'salmon'
kallisto_folder = "kallisto"
bowtie_bin = "bowtie2" #bowtie2 binary
#----
blast_folder = '/home/phyllobates/NCBI_Tools/ncbi-blast-2.7.1+/'
custom_blast_folder = '/home/phyllobates/annotation_libraries/Uniprot_library' #custom blast library location
blast_name = "uniprot_db" #name of custom blast library; name in db folder
evalue = "1e-10"
trim_length = 300

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