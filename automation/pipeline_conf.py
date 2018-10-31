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

#name of custom blast library
blast_name = "name_of_library"

#final file path for assembly output
filename = "final.contigs.fa"

#location of forward and reverse reads
forwreads = '/home/litoria/Desktop/media/genomes/Corytophanes_percarinatus_unpublished/hemipenes/3M-Nathan_S3_R1_001.fastq.gz'
revreads = '/home/litoria/Desktop/media/genomes/Corytophanes_percarinatus_unpublished/hemipenes/3M-Nathan_S3_R2_001.fastq.gz'
