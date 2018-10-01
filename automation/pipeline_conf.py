#input file location
input_srx = "path/to/file"

#custom blast library location
custom_location = "/home/litoria/Assembly_Tools/Uniprot_library" 	#this might need to be removed if above input works

#name of custom blast library
blast_name = "name_of_library"

#final file path for assembly output
filename = "final.contigs.fa"

#location of megahit bin (if in $PATH just place the name of the bin)
megahit = "megahit"
blast_bin = "/home/litoria/NCBI_Tools/ncbi-blast-2.7.1+/bin/blastp"

assembly_dir = f"{srx}.megahit_asm"

working_dir = "/home/litoria/Documents/test"
