set -x
input=$(yad --title="Megahit" --text="Select input files:" \
	--form \
	--field="Working Directory":DIR \
	--field="Forward Reads":FL \
	--field="Reverse Reads":FL)

wdinput="$(echo "$input"| cut -d '|' -f 1)"
forwreads="$(echo "$input"| cut -d '|' -f 2)"
revreads="$(echo "$input" | cut -d '|' -f 3)"

#extract filenames w/ extension
forwinput=$( echo "${forwreads##*/}")
revinput=$( echo "${revreads##*/}")

#extract filenames w/o extensions
forwfilename=$(echo "$forwinput" | cut -d '.' -f 1)
revfilename=$(echo "$revinput" | cut -d '.' -f 1)

#run megahit in appropriate mode
cd "$wdinput"
megahitfolder=$".megahit_asm" 
count=$(ls -R | grep "\.fastq$" | wc -l)
count2=$(ls -R | grep "\.fastq.gz$" | wc -l)
	if [ "$count" -eq 2 ] || [ "$count2" -eq 2 ]; then
	(megahit -1 "$forwinput" -2 "$revinput" -o "$PWD$megahitfolder")
	mv -v "$PWD$megahitfolder" "$wdinput"
	else
		if [ "$count" == 1 ]; then
		(megahit -r "$forwinput" -o "$PWD$megahitfolder")
		mv -v "$PWD$megahitfolder" "$wdinput"
		else
			if [ "$count" == 0 ] || [ "$count" -gt 2 ]; then
			(zenity --error --title="An Error Occurred" --text="no FASTQ files were found")
			fi
		fi
	fi
