#!/usr/bin/env python3
import sys
sys.path.insert(0, '/home/litoria/Documents/git/transcriptome_assembly/automation/')

from automation_script import *

forwreads = []
revreads = []

#change line below to pick which modules to run
def main():
    manual_input()
    # trimmomatic()
    # megahit1()
    # abyss()
    # spades()
    # transrate()
    # annotation()
    # pull_matches()
    # pull_matches_fast()
    rsem()

if __name__ == '__main__':
    main()


#input file location
input_srx = "path/to/file"

#binary locations
megahit_bin = "megahit"
abyss_bin = "abyss-pe"
spades_bin = "/home/litoria/Assembly_Tools/SPAdes-3.12.0/bin/spades.py" #location of spades binary on system
transrate_bin = "/home/litoria/Assembly_Tools/transrate-1.0.3-linux-x86_64/transrate"
blast_bin = "/home/litoria/NCBI_Tools/ncbi-blast-2.7.1+/bin/blastn"
trimmomatic_jar = "/home/litoria/Assembly_Tools/Trimmomatic-0.38/trimmomatic-0.38.jar"
rsem_loc = '/home/litoria/Assembly_Tools/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl'
bowtie_bin = "bowtie2"

#custom blast library location
custom_location = "/home/litoria/Assembly_Tools/blastn_mitochondria"

#name of custom blast library; name in db folder
blast_name = "frog_mito_db"
