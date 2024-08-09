import pandas as pd
import boto3
import os

# Raw reads folder path relative to experiment folder
raw_path = "seq/00_raw_reads"

# Create a directory to store sequencing data
if not os.path.isdir(raw_path):
    os.makedirs(raw_path)

# Load metadata
meta = pd.read_csv("metadata/metadata.csv", index_col = "seqid")

# Global S3/boto3 info & connection
s3 = boto3.client('s3')
BUCKET_NAME = "guillaume-urtecho"

# Downloads reads for this experiment into raw_path from YH S3 master folder
# Used as an apply function on the metadata df
def download_reads(row):
    # Constructing local paths to dl seq files to
    seqid = row.name
    R1_path = "{}/{}_R1.fastq.gz".format(raw_path, seqid)
    R2_path = "{}/{}_R2.fastq.gz".format(raw_path, seqid)

    # Constructing S3 path to seq files in GU folder
    S3_R1_path = row.S3_R1_path
    S3_R2_path = row.S3_R2_path

    # Downloading from S3 to local folder 
    s3.download_file(BUCKET_NAME, S3_R1_path, R1_path)
    s3.download_file(BUCKET_NAME, S3_R2_path, R2_path)

# Start dl
meta.apply(download_reads, axis=1)