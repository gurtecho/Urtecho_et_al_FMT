
###Prepare directory###
-Make all the experiment-seq directories & put raw data with correct naming e.g. {seqid}_R1.fastq.gz

If need to rewrite files names from basespace, can use this loop:

for file in *.fastq.gz; do
    new_name=$(echo "$file" | sed 's/_S[0-9]\+_L001_/_/; s/_001//')
    mv "$file" "$new_name"
done


#Load up thomas's AMI on T2x2.large
Tum3, use personal key, set traffic to allow current IP address (if not coding on the move)


##Log into AWS

ssh -i ~/Google_Drive/WangLab/guillaume.pem ec2-user@<IPV4 address here>

#Synch data from s3 bucket
#example:  aws s3 sync s3://guillaume-urtecho/20230613_MapsBulk_ec2/ --endpoint-url https://s3.us-east-1.amazonaws.com .

aws s3 sync <S3 URL> --endpoint-url https://s3.us-east-1.amazonaws.com .


#Install Snakemake
conda install -c bioconda -c conda-forge snakemake

#run snakemake in directory with snakemake file, this was with t2x2.large
snakemake --cores 4


