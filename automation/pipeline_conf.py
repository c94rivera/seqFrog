#!/usr/bin/env python3
import sys
sys.path.insert(0, '/home/phyllobates/Documents/git/transcriptome_assembly/automation')
from automation_script import *
forwreads = []
revreads = []

#binary locations
megahit_bin = "megahit" #megahit binary
abyss_bin = '/home/phyllobates/Assembly_Tools/abyss-2.1.5/ABYSS/abyss.cc' #abyss binary
spades_bin = '/home/phyllobates/Assembly_Tools/SPAdes-3.13.0-Linux/bin/spades.py' #spades binary
transrate_bin = '/home/phyllobates/Assembly_Tools/transrate-1.0.3-linux-x86_64/transrate' #transrate binary
trimmomatic_jar = '/home/phyllobates/Assembly_Tools/Trimmomatic-0.38/trimmomatic-0.38.jar' #trimmomatic binary
#rsem_loc should be the location of the "align_and_estimate_abundance.pl" script within your Trinity install
rsem_loc = '/home/litoria/Assembly_Tools/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl' #rsem location
salmon_bin = 'salmon'
bowtie_bin = "bowtie2" #bowtie2 binary
#----
blast_bin = '/home/phyllobates/NCBI_Tools/ncbi-blast-2.7.1+/bin/blastp' #blast binary; use blastn for nucleotide and blastp for protein annotation
custom_location = '/home/phyllobates/annotation_libraries/Uniprot_library' #custom blast library location
blast_name = "uniprot_db" #name of custom blast library; name in db folder

#change line below to pick which modules to run
def main():
    manual_input()
    trimmomatic()
    megahit()
    abyss()
    spades()
    transrate()
    annotation()
    blastx()
    pull_matches() #use only if pull_matches_fast doesn't work
    pull_matches_fast()
    rsem()

if __name__ == '__main__':
    main()


#input file location
input_srx = "path/to/file"
