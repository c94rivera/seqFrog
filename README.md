# seqFrog

A modular pipeline for transcriptomic work in non-model organisms.

If using this pipeline please properly cite this git repository:

`Rivera, Christopher, seqFrog, (2019), GitHub repository, https://github.com/c94rivera/seqFrog`


### Get Started
The following programs need to be installed in your $PATH variable by following their full installation instructions. This should be done FIRST with the installation of the pipeline.
Run the following commands then proceed with Samtools and Bowtie2 installation instructions
```
sudo apt install gcc
sudo apt install make
sudo apt install libbz2-dev
sudo apt install libncurses5-dev
sudo apt install zlib1g-dev
sudo apt install libncursesw5-dev
sudo apt install liblzma-dev
sudo apt install libcurl4-openssl-dev
sudo apt install libssl-dev
```

* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools](http://www.htslib.org/)

Download the repository to your computer. The following programs must be installed and their folder PATH setup in seqFrog_conf.py file. More details can be found within Methods.

* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Megahit](https://github.com/voutcn/megahit)
* [ABySS](https://github.com/bcgsc/abyss)
* [SPAdes](http://cab.spbu.ru/software/spades/)
* [Transrate](http://hibberdlab.com/transrate/)
* [BLAST+ command line applications](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
* [RSEM](https://deweylab.github.io/RSEM/)
* [Salmon](https://combine-lab.github.io/salmon/)
* [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
* [Kallisto](https://pachterlab.github.io/kallisto/)


The following Python3 modules are necessary for the use of seqFrog:

pandas `pip3 install pandas`
Biopython `pip3 install biopython`
tqdm `pip3 install tqdm`

You will also need a working version of Java

---
# Follow Along Guide
### Prep-work
1. Direct to all necessary program folders
	* Open `seqFrog_conf.py` in text-editor of choice
	* Put the full length path to all the required program folders next to their similarly named variables, in quotes `lines 11-25` (relative paths will not work!)
		* if your programs are installed into system or $PATH variables you can just use the name used to call it in the terminal in quotes
	* Put the full length path to the custom blast library being used (relative path will not work!)

```python3
megahit_folder = "/megahit" #megahit binary
abyss_folder = '/home/phyllobates/Assembly_Tools/abyss-2.1.5' #abyss binary
spades_folder = '/home/phyllobates/Assembly_Tools/SPAdes-3.13.0-Linux' #spades binary
```


2. Define steps that should be taken within pipeline
	* The functions after `def main():` can be commented out to exclude steps from the pipeline
		* the pipeline is designed to be run from start to finish, however, individual steps can be run in some instances

```python3
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
```

### seqFrog Initiation
1. Open the folder you wish to work within (usually where the reads are located)
2. Open a terminal window within this folder
3. Call the `seqFrog_conf.py` file
	this can be done by typing out the path to the file or dragging and dropping the file into the terminal
4. Fill in required flags (detailed below)
5. Fill in optional flags if desired or running only a portion of seqFrog
6. Wait for output!

#### Flags
	Required input:
		"-f" or "--forward-reads" = Forward Reads
		"-r" or "--reverse-reads" = Reverse Reads

	Optional input:
		"-c" or "--contigs" = Contig file (used if skipping assemblies/assemblies already done)
		"-b" or "-blast-file" = Blast file for comparison (used if skipping annotation within seqFrog/annotation already done)
		"-sp" or "--species" = Species name for organism being analysed (used for naming of final output files)
		"-t" or "--tissue" = Tissue type for organism being analysed (used for naming of final output files)

	Do not use special symbols in file names or expression analysis will not function properly.

#### Example
```bash
seqFrog_conf.py -f fwdreads.fastq -r revreads.fastq -sp Dendrobates_pumilio -t brain
```
---
### How to make a custom blast library
coming soon...
