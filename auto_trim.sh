#ask user for inputs
input=$(yad --title="Trimmomatic" --text="Select input files:" \
	--form \
	--field="Working Directory":DIR \
	--field="Forward Reads":FL \
	--field="Reverse Reads":FL)

echo $input

# wdinput=$(cut -d '\n' -f1 /tmp/entries)
# forwreads=$(cut -d '|' -f2 /tmp/entries)
# revreads=$(cut -d '|' -f3 /tmp/entries)

#	java -jar trimmomatic-0.38.jar PE -threads 10 -phred33 $forwinput $revinput "$filename1"_paired.fastq "$filename1"_unpaired.fastq "$filename2"_paired.fastq "$filename2"_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
