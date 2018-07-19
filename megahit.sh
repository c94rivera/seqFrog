###this is to run it in verbose mode for debugging
set -x

#ask user for inputs
input=$(yad --title="Megahit" --text="Select input files:" \
	--form \
	--field="Bulk File":FL \
	--field="Working Directory":DIR \
	--field="Forward Reads":FL \
	--field="Reverse Reads":FL )

#parse $input and seperate different inputs
bulkfile="$(echo "$input"| cut -d '|' -f 1)"
wdinput="$(echo "$input"| cut -d '|' -f 2)"
forwreads="$(echo "$input"| cut -d '|' -f 3)"
revreads="$(echo "$input" | cut -d '|' -f 4)"

#if there is a bulk file then run megahit automatically
if [ -z "$bulkfile" ]; then
	xargs mkdir "$bulkfile"


		
	else
		cd $newfolder
		datapath_f=$(zenity --file-selection --multiple --title="Select Forward Reads")
		datapath_r=$(zenity --file-selection --multiple --title="Select Reverse Reads")

fi

#decides which programs to run based on user input
userinput2=$(zenity --list --checklist --title="Options"\
    --text="Select your features"\
    --column="Use"\
    --column="Feature"\
    TRUE "Megahit Assembly"\
    TRUE "BLAST Annotation"\
    FALSE "Custom BLAST Library")
megahitdecision="$(echo "$userinput2"| cut -d '|' -f 1)"
blastdecision="$(echo "$userinput2"| cut -d '|' -f 2)"
blastlibrarydecision="$(echo "$userinput2" | cut -d '|' -f 3)"
#
#
srx1=$"_1.fastq" 
srx2=$"_2.fastq" 
for d in ./*/ ; do (cd "$d" &&
foldername=$"${PWD##*/}" 
srxnum=$( echo "${PWD##*/}" | awk -F'[_]' '{print $3}') 
cd .. 
if [ "$megahitdecision" == "Megahit Assembly" ]; then
	megahitfolder=$".megahit_asm" 
	count=$(ls -R |grep "\.fastq$" |wc -l)
	count2=$(ls -R |grep "\.fastq.gz$" |wc -l)
		if [ "$count" -eq 2 ] || [ "$count2" -eq 2 ]; then
		(megahit -1 "$datapath_f" -2 "$datapath_r" -o "$PWD$megahitfolder")
		else
			if [ $count == 1 ]; then
			(megahit -r "$datapath_f" -o "$PWD$megahitfolder")
			else
				if [ $count == 0 ]; then
				zenity –error –title=”An Error Occurred” –text=”no FASTQ files were found”
				fi
			fi
		fi
	mv -v "$PWD$megahitfolder" $foldername 
fi
###add interface for selecting blast library or using standard one
###add automatic blast annotation if library is not selected
if [ "$blastlibrarydecision" == "Custom BLAST Library" ]; then
	cp -R $blastlibrarypath . 
	cd "${PWD##*/}".megahit_asm 
	cp final.contigs.fa .. 
	cd .. 
	mv final.contigs.fa blastn_mitochondria 
	cd blastn_mitochondria 
	blastn -query final.contigs.fa -db frog_mito_db -evalue 0.01 -outfmt 7 -out "$foldername"_megahit.table 
	cd ..
	mkdir $foldername 
	cd blastn_mitochondria 
	cp final.contigs.fa ../"$foldername" 
	cp "$foldername"_megahit.table ../"$foldername"
fi);
done