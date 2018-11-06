#input file location
input_srx = "path/to/file"

#binary locations
megahit_bin = "megahit"
abyss_bin = "abyss-pe"
spades_bin = "/home/litoria/Assembly_Tools/SPAdes-3.12.0/bin/spades.py" #location of spades binary on system
transrate_bin = "/home/litoria/Assembly_Tools/transrate-1.0.3-linux-x86_64/transrate"
blast_bin = "/home/litoria/NCBI_Tools/ncbi-blast-2.7.1+/bin/blastp"

#custom blast library location
custom_location = "/home/litoria/Assembly_Tools/Uniprot_library" 	#this might need to be removed if above input works

#name of custom blast library; name in db folder
blast_name = "uniprot_db"

#contig file for blast annotation, only use is not running assembly first
# contig_file = []

#location of forward and reverse reads
