#!/usr/bin/env python3
import sys
sys.path.insert(0, '/home/phyllobates/Documents/git/transcriptome_assembly/automation')

from automation_script import *

forwreads = []
revreads = []

#change line below to pick which modules to run
def main():
    manual_input()
    trimmomatic()
    megahit1()
    abyss()
    spades()
    transrate()
    annotation()
    # pull_matches() #use only if pull_matches_fast doesn't work
    pull_matches_fast()
    rsem()

if __name__ == '__main__':
    main()


#input file location
input_srx = "path/to/file"

#binary locations
megahit_bin = "megahit"
abyss_bin = '/home/phyllobates/Assembly_Tools/abyss-2.1.5/ABYSS/abyss.cc'
spades_bin = '/home/phyllobates/Assembly_Tools/SPAdes-3.13.0-Linux/bin/spades.py' #location of spades binary on system
transrate_bin = '/home/phyllobates/Assembly_Tools/transrate-1.0.3-linux-x86_64/transrate'
blast_bin = '/home/phyllobates/NCBI_Tools/ncbi-blast-2.7.1+/bin/blastp'
trimmomatic_jar = '/home/phyllobates/Assembly_Tools/Trimmomatic-0.38/trimmomatic-0.38.jar'
#rsem_loc should be the location of the "align_and_estimate_abundance.pl" script within your Trinity install
rsem_loc = '/home/phyllobates/Assembly_Tools/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl'
bowtie_bin = "bowtie2"

#custom blast library location
custom_location = '/home/phyllobates/annotation_libraries/Uniprot_library'

#name of custom blast library; name in db folder
blast_name = "uniprot_db"
