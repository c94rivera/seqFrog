set -x
wd=$(zenity --file-selection --directory --title="Select Working Directory")
cd $wd
datapath_f=$(zenity --file-selection --multiple --title="Select Forward Reads")
datapath_r=$(zenity --file-selection --multiple --title="Select Reverse Reads")

megahitfolder=$".megahit_asm" 
count=$(ls -R |grep "\.fastq$" |wc -l)
count2=$(ls -R |grep "\.fastq.gz$" |wc -l)
	if [ "$count" -eq 2 ] || [ "$count2" -eq 2 ]; then
	(megahit -1 "$datapath_f" -2 "$datapath_r" -o "$PWD$megahitfolder")
	else
		if [ "$count" == 1 ]; then
		(megahit -r "$datapath_f" -o "$PWD$megahitfolder")
		else
			if [ "$count" == 0 ] || [ "$count" -gt 2 ]; then
			(zenity --error --title="An Error Occurred" --text="no FASTQ files were found")
			fi
		fi
	fi
mv -v "$PWD$megahitfolder" $foldername