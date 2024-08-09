

#Load up tum AMI on aws

#Install anvio

conda update conda
conda create -y --name anvio-8 python=3.10
conda activate anvio-8

#install mamba for fast dependencies
conda install -y -c conda-forge mamba


mamba install -y -c conda-forge -c bioconda python=3.10 \
        sqlite prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
        blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
        trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
        r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript


# try this, if it doesn't install, don't worry (it is sad, but OK):
#mamba install -y -c bioconda fastani

curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz \
        --output anvio-8.tar.gz

#maybe mamba can work
conda install -c conda-forge -c bioconda gcc
pip install anvio-8.tar.gz

#######Colonizer vs non-colonizer analysis##########

# download donor A data pack:
curl -L https://ndownloader.figshare.com/files/27452192 \
     -o FMT_DONOR_A_AND_RECIPIENTS.tar.gz


tar -zxvf FMT_DONOR_A_AND_RECIPIENTS.tar.gz

cd FMT_DONOR_A_AND_RECIPIENTS

# download additional input files
for file in scg-cov-DA.txt \
            scg-cov-DB.txt \
            subpop-comp-DA.txt \
            subpop-comp-DB.txt \
            subpop-num-DA.txt \
            subpop-num-DB.txt \
            metadata-donor.txt \
            metadata-recipient.txt \
            metadata-transplants.txt \
            detection-FMT-DA.txt \
            detection-FMT-DB.txt \
            detection-global-by-country-DA.txt \
            detection-global-by-country-DB.txt;
do curl -L https://merenlab.org/data/fmt-gut-colonization/files/${file} \
        -o ${file};
done


# download the script
curl -L https://merenlab.org/data/fmt-gut-colonization/files/determine-colonization.py \
     -o determine-colonization.py


python determine-colonization.py


#Do CAZyme analysis, default uses CAZydb v11 so do it with version 12; gets you contigs_cazymes.txt
anvi-setup-cazymes --cazyme-version V12
anvi-migrate CONTIGS.db --migrate-safely
anvi-run-cazymes -c CONTIGS.db -T 30
anvi-export-functions -c CONTIGS.db --annotation-sources CAZyme -o contigs_cazymes.txt

#Export binning information, gets you my_bins.txt
anvi-migrate PROFILE.db --migrate-safely
aanvi-summarize -c CONTIGS.db -p PROFILE.db -o SUMMARY_DIR --init-gene-coverages -T 30


#Export bin calls, gets you 
anvi-export-gene-calls -c CONTIGS.db --skip-sequence-reporting -o gene_calls.txt --gene-caller prodigal



#Run this python script to merge

```
import pandas as pd

# Load the datasets
genes_in_splits = pd.read_csv("../genes_in_splits.tsv", sep="\t")
contigs_cazymes = pd.read_csv("../contigs_cazymes.txt", sep="\t")
my_bins = pd.read_csv("../my_bins.txt", sep="\t", names=["split", "MAG"])

# Merge contigs_cazymes with genes_in_splits
merged_df = pd.merge(contigs_cazymes, genes_in_splits, on="gene_callers_id", how="left")

# Merge the result with my_bins
final_df = pd.merge(merged_df, my_bins, on="split", how="left")

# Print out the number of unique splits in both genes_in_splits and my_bins
print(f"Unique splits in genes_in_splits: {genes_in_splits['split'].nunique()}")
print(f"Unique splits in my_bins: {my_bins['split'].nunique()}")

# Find splits in genes_in_splits that are not in my_bins
missing_splits = set(genes_in_splits['split']) - set(my_bins['split'])
print(f"Number of splits in genes_in_splits not in my_bins: {len(missing_splits)}")

# Optionally, print some of these missing splits to inspect
print("Some missing splits:", list(missing_splits)[:10])

#drop NA in MAG
final_df.dropna(subset=['MAG'], inplace=True)

# Print the first few rows to check the result
print(final_df.head())

# Save the final dataframe to a CSV file
final_df.to_csv("combined_cazymes_MAGs.csv", index=False)

```







