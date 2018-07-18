function select_files_forward 
{ 
    yad --title "Select forward read files" --file-selection --multiple
}
export -f select_files_forward
function select_files_reverse 
{ 
    yad --title "Select reverse read files" --file-selection --multiple
}
export -f select_files_reverse
function select_working_directory 
{ 
    yad --title "Select working directory" --file --directory
}
export -f select_working_directory


#	java -jar trimmomatic-0.38.jar PE -threads 10 -phred33 $forwinput $revinput "$filename1"_paired.fastq "$filename1"_unpaired.fastq "$filename2"_paired.fastq "$filename2"_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

input=$(yad --form --width=250 --text="Trimmomatic" \
	--field="Select Working Directory":btn "bash -c select_working_directory" \
	--field="Select Forward Reads":btn "bash -c select_files_forward" \
	--field="Select Reverse Reads":btn "bash -c select_files_reverse")
#seperate the variables into individual variables
wdinput="$(echo "$input"| cut -d ' ' -f 1)"
forwreads="$(echo "$input"| cut -d ' ' -f 2)"
revreads="$(echo "$input"| cut -d ' ' -f 3)"

echo $wdinput
echo $forwreads
echo $revreads

#select working directory
#get name of file for output names
#yad --form --title="Trimmomatic" \
#	--button="Select Working Directory:bash -c select_working_directory" \
#	--button="Select Forward Reads:bash -c select_files_forward" \
#	--button="Select Reverse Reads:bash -c select_files_reverse"