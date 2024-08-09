RUN=$1

#This will create a folder named 'fastq' and download all reads to it. 
#You can alter the output directory name by altering the name after the -o flag in the bs download command

###  USAGE ###
#activate script with:
# chmod +x DownloadBasespace.sh
#
#run command with:
#
#./DownloadBasespace.sh <run ID #> 
#
#You can find run ID # through basespace by going to your recent projects and copying the number embedded in the from the URL. See example url:
#Project URL: https://basespace.illumina.com/projects/396032646/about
#Example run command: ./DownloadBasespace.sh 396032646

												###########  Setup Basespace ###########

echo -e "Do you need to install basespace CLI? Plz respond (yes/no) \n"
read response

# Check the user's response and execute different sections of code
if [ "$response" == "yes" ]; then
    echo 'Installing Basespace CLI, get ready to log in! :)'

#Get BaseSpace CLI tool
wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O $HOME/bin/bs

#Enable read/write permissions
chmod u+x $HOME/bin/bs

elif [ "$response" == "no" ]; then
    echo "No need to install! :)"
else
    echo "Invalid response. Please enter 'yes' or 'no'."
    exit 1  # Exit the script with a non-zero status code
fi

#Authenticate account 
echo -e "Follow prompts to log into basespace account \n"
bs auth




											###########  Downloading files  ###########

echo -e "Downloading fastq files \n"
mkdir -p fastq
bs download project -i $RUN -o fastq --extension=fastq.gz





											###########  Extract from ds directories  ###########

echo 'Extracting sequencing files from basespace \n'

cd fastq


for i in *ds.* 
do	
	cd $i
	mv *.fastq.gz ..
	cd ..
done

rm -r *ds.*



											 ###########  Merge Lanes  ###########

echo -e "Do you want to combine sample files across lanes? Plz respond (yes/no) \n"
read response

# Check the user's response and execute different sections of code
if [ "$response" == "yes" ]; then
    echo "Combining lanes :)"

			for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
			    do echo "Merging R1"
					echo $i
						cat "$i"*_R1.fastq.gz > "$i"_ME_R1.fastq.gz
			       	echo "Merging R2"
						cat "$i"*_R2.fastq.gz > "$i"_ME_R2.fastq.gz
			done;
elif [ "$response" == "no" ]; then
    echo "Not combining, good luck! :)"
else
    echo "Invalid response. Please enter 'yes' or 'no'."
        exit 1  # Exit the script with a non-zero status code
fi


											 ###########  Fix names  ###########
echo 'Fixing file names \n'

for file in *_S*R1_001.fastq.gz; do mv "$file" "${file/_S*_001.fastq.gz/}_R1.fastq.gz"; done
for file in *_S*R2_001.fastq.gz; do mv "$file" "${file/_S*_001.fastq.gz/}_R2.fastq.gz"; done

# Continue with the rest of your script
echo "Script complete. Have a nice day!  \n"


echo -e " HHH  HHH  IIII"
echo -e " HHH  HHH   II "
echo -e " HHHHHHHH   II "
echo -e " HHH  HHH   II "
echo -e " HHH  HHH  IIII\n\n"
echo -e "WW    WW   AAA   NNN   NN  GGGGGG  \nWW    WW  AAAAA  NNNN  NN GGG      \nWW WW WW AA   AA NN NN NN GGG  GGG \nWW WW WW AAAAAAA NN  NNNN GGG   GG \n WWWWWW  AA   AA NN   NNN  GGGGGG  "
echo -e "\n"
echo -e " LL       AAA   BBBB   \n LL      AAAAA  BB   B  \n LL     AA   AA BBBBBB  \n LL     AAAAAAA BB   BB \n LLLLLL AA   AA BBBBBB  \n\n"


