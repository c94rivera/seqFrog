#ask user for inputs
input=$(yad --title="Trimmomatic" --text="Select input files:" \
	--form \
	--field="Working Directory":DIR \
	--field="Forward Reads":FL \
	--field="Reverse Reads":FL \
	--field="Trimmomatic jar file":FL)

#parse $input and seperate different inputs
wdinput="$(echo "$input"| cut -d '|' -f 1)"
forwreads="$(echo "$input"| cut -d '|' -f 2)"
revreads="$(echo "$input" | cut -d '|' -f 3)"
trimmo="$(echo "$input" | cut -d '|' -f 4)"

#extract filenames w/ extension
forwinput=$( echo "${forwreads##*/}")
revinput=$( echo "${revreads##*/}")

#extract filenames w/o extensions
forwfilename=$(echo "$forwinput" | cut -d '.' -f 1)
revfilename=$(echo "$revinput" | cut -d '.' -f 1)

#run trimmomatic with default settings and input files from above
cd $wdinput
java -jar "$trimmo" \
	-threads 10 -phred33 \
	$forwinput \
	$revinput \
	"$forwfilename"_paired.fastq \
	"$forwfilename"_unpaired.fastq \
	"$revfilename"_paired.fastq \
	"$revfilename"_unpaired.fastq I \
	LLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
