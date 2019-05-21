# seqFrog

A modular pipeline for transcriptomic work in non-model organisms.

If using this pipeline please properly cite this git repository:

`Rivera, Christopher, seqFrog, (2019), GitHub repository, https://github.com/c94rivera/seqFrog`


### Get Started
Download the repository to your computer by click 'clone or download' button on the top right
The following dependecies need to be installed on your computer. This should be done FIRST, before the installation of the pipeline.
```
sudo apt install gcc
sudo apt install make
sudo apt install cmake
sudo apt install libbz2-dev
sudo apt install libncurses5-dev
sudo apt install zlib1g-dev
sudo apt install libncursesw5-dev
sudo apt install liblzma-dev
sudo apt install libcurl4-openssl-dev
sudo apt install libssl-dev
sudo apt install libsparsehash-dev
```
Install bowite2 and samtools using the following commands (recommended) or the instructions on their website.
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) `sudo apt install bowtie2`
* [Samtools](http://www.htslib.org/) `sudo apt install samtools`

The following program binary files must be installed and their folder PATH setup in seqFrog_conf.py file.
It is recommended that all these programs be downloaded and extracted into the same folder for easy setup and maintenance.

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)  
Download binary for version 0.39 (not source code)

[Megahit](https://github.com/voutcn/megahit)  
Click on releases and download latest stable version (not source code)

[ABySS](https://github.com/bcgsc/abyss)  
Click releases and download latest stable version 2.1.5 (not source code)  
Right click in the ABySS folder and open a new terminal. Do the following commands in the terminal to install ABySS  
```
./configure
make
sudo make install
```

[SPAdes](http://cab.spbu.ru/software/spades/)  
Download binary for linux version 3.13.0

[Transrate](http://hibberdlab.com/transrate/)  
Click on installation and download linux version 1.0.3

[BLAST+ command line applications](https://www.ncbi.nlm.nih.gov/books/NBK279671/)  
Click on the FTP link, click LATEST and download the XXX_linux.tar.gz

[RSEM](https://deweylab.github.io/RSEM/)  
Click 'latest version' to download the source code
Right click in the RSEM folder and open a new terminal. Do the following commands in the terminal to install RSEM
```
make
sudo make install
```

[Salmon](https://combine-lab.github.io/salmon/)  
Click on binaries, on the github page download the binaries for version 0.13.1

[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)  
Follow links for download and download the source code for version 2.8.4
Right click in the Trinity-RNAseq folder and open a new terminal. Do the following commands in the terminal to install Trinity
```
make
sudo make install
```

* [Kallisto](https://pachterlab.github.io/kallisto/)  
Click on downloads link, scroll to releases and download Linux version 0.45.0

The following Python3 modules are necessary for the use of seqFrog:

pandas `pip3 install pandas`  
Biopython `pip3 install biopython`  
tqdm `pip3 install tqdm`  
Remove any older versions of these modules

You will also need a working version of Java

### It is recommended that you restart the computer after installing all the programs.

---
# Follow Along Guide
### Prep-work
1. Direct to all necessary program folders
	* Open `seqFrog_conf.py` in text-editor of choice
	* Put the full length path to all the required program folders next to their similarly named variables, in quotes `lines 11-25` (relative paths will not work!)
		* The easiest way to get the full length path on linux is to drag the folder into a terminal window and the full length path will be inserted into the terminal window. Copy and past this into the appropriate part of `seqFrog_conf.py`
		* if bowtie2 was installed as root then use 'bowtie2', otherwise, direct to the folder.
	* Put the full length path to the custom blast library being used

```python3
megahit_folder = '' #megahit binary
abyss_folder = '/home/phyllobates/Assembly_Tools/abyss-2.1.5' #abyss binary
spades_folder = '/home/phyllobates/Assembly_Tools/SPAdes-3.13.0-Linux' #spades binary
```


2. Define steps that should be taken within pipeline
	* The functions after `def main():` can be commented out () to exclude steps from the pipeline
		* the pipeline is designed to be run from start to finish, however, individual steps can be run in some instances

### Function Descriptions
	manual_input() -> takes the input from the terminal, CANNOT BE TURNED OFF
	trimmomatic() -> runs Trimmomatic to clean reads
	megahit() -> runs Megahit for de novo reconstruction
	abyss() -> runs ABySS for de novo reconstruction
	spades() -> runs SPAdes for de novo reconstruction
	transrate() -> runs Transrate to combine reconstructions and give statistics
	trim_contigs() -> removes contigs below a certain base pair length
	blastn() -> runs BlastN search against custom library
	blastp() -> runs BlastP search against custom library
	blastx() -> runs BlastX search against custom library
	pull_matches() -> annotates reconstructions with Blast data
	pull_matches_fast() -> annotates reconstructions with Blast data using all available cores on computer
	rsem() -> runs RSEM for expression analysis
	salmon() -> runs Salmon for expression analysis
	kallisto() -> runs Kallisto for expression analysis

For example the following code skips blastn(), blastp(), pull_matches().
```python3
def main():
    manual_input()
    trimmomatic()
    megahit()
    abyss()
    spades()
    transrate()
    trim_contigs()
    # blastn()
    # blastp()
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
