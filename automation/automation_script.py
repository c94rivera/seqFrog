import os
import shutil
import glob
import subprocess

#this script assumes that all the necessary bin files are in PATH
#it also is designed for megahit as the only assembler

########## User input
custom_location = "/home/litoria/Assembly_Tools/Uniprot_library" 	#location of custom blast library
blast_name = "uniprot_db"		#name of custom blast library
filename = "final.contigs.fa"	#name of output from assembler
input_srx = "/home/litoria/Documents/test/foldernames.txt"	#name of file with all SRX numbers and sample names
blast_bin = "/home/litoria/NCBI_Tools/ncbi-blast-2.7.1+/bin/blastp"
working_dir = "/home/litoria/Documents/test"
assembly_dir = f"{srx}.megahit_asm"
output_file = f"{subdir}_megahit.table"



########## Functions


#make folders based on strings in input file
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


#run fastqdump on specific srx number
def fastqdump(srx):
	subprocess.call(["fastq-dump", "--dumpbase", "--defline-seq", "'@$sn[_$rn]/$ri'", "--split-files", srx])



#run megahit in single or paired end mode depending on amount of input files found
def megahit(srx):
		if len(glob.glob("*.megahit_asm")) >= 1:
			print("Megahit assembly folder already present")
		elif len(glob.glob("*.fastq")) == 2:
			subprocess.call(["megahit", "-1", f"{srx}_1.fastq", "-2", f"{srx}_2.fastq", "-o", f"{srx}.megahit_asm"])
		elif len(glob.glob("*.fastq")) == 1:
			subprocess.call(["megahit", "-r", f"{srx}_1.fastq", "-o", f"{srx}.megahit_asm"])
		else:
			print("Error: Input files were not found")


#copy custom blast library into working directory
# def customblast():
# 	try:
# 		if not os.path.exists(os.getcwd() + "/blast_annotation"):
# 			os.mkdir(os.getcwd() + "/blast_annotation")
# 			shutil.copytree(custom_location, os.getcwd() + "/blast_annotation")
# 	except OSError:
# 		print("Blast directory could not be copied.")

#copy final.contigs.fa into customblast folder and run annotation
def annotation():
	contig_file = None
	blast_code = [f"{blast_bin}", "-query", f"{contig_file}", "-db", f"{blast_name}", "-evalue", "0.01", "-outfmt", "7", "-out", f"{subdir}_megahit.table"]
	for dirs in os.listdir():
		if assembly_dir in dirs:
			print(f"Assembly folder found: running annotation with {blast_name} library")
			contig_file = (os.getcwd() + "/" + assembly_dir + "/" + filename)
			current_dir = os.getcwd()
			os.chdir(custom_location)
			print(blast_code)
			subprocess.call([f"{blast_bin}", "-query", f"{contig_file}", "-db", f"{blast_name}", "-evalue", "0.01", "-outfmt", "7", "-out", f"{subdir}_megahit.table"])
			os.copy()
			os.chdir(current_dir)
			#add something that will tell the user it is still running every few seconds
		else:
			print("Error: Assembly folder not found")


#create final output folder and move output files into this folder
def finalfolder():
	output_file = f"{subdir}_megahit.table"
	os.mkdir(subdir)
	os.chdir(subdir)
	shutil.copy(custom_location + "/" + filename, ".")
	shutil.copy(os.getcwd() + "/" + assembly_dir + filename, "." + subdir)


########## End of Functions


#make folders based on the lines within a textfile
os.chdir(working_dir)
makefolders(input_srx)

subdirs = next(os.walk("."))[1]
print(subdirs)												#remove after testing
for subdir in subdirs:
	#grab srx from folder name
	print(subdir)
	srx = subdir.split("_")[2]
	print(srx)
	#open folder
	os.chdir(os.getcwd() + "/" + subdir)
	print(os.getcwd())										#remove after testing
	#run fastqdump
	fastqdump(srx)
	#run megahit
	megahit(srx)
	#prepare folders for blast blast_annotation
	# customblast()
	#run annotation on custom blast library
	annotation()
	#make additional folder and copy final files into it
	finalfolder()
	#leave current folder to continue loop
	os.chdir("..")
	os.chdir("..")
