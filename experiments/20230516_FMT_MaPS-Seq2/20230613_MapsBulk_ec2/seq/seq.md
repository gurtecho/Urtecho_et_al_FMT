##
## 16S Pipeloine

### 2021-05-30
#### Running Full 16S Pipeline across all Experiments
* Running 16S pipeline across experiments:
    * `../experiments/2017-10-10_pairwise_cohousing`
    * `../experiments/2018-09-30_JxE_FMT`
    * `../experiments/2018-11-13_JxE_FMT02`
    * `../experiments/2018-11-14_JxE_BCT01`
    * `../experiments/2019-07-11_JxE_BCT02`
    * `../experiments/2021-03-18_JxE_FMT03`
    * `../experiments/2021-05-12_envigo_culture_pilot`
* Using snakemake
* Constructing pipeline as I go
* Ran out of memory for Derep! Will try again tmrw w/ larger EC2 instance



### 2021-06-31
#### Reconfiguring EC2 Instance
* Launched an c5.12xlarge instance
* Connected : 
* Pulled down data
    * `mkdir fmt`
    * `cd fmt`
    * `mkdir experiments; mkdir seq`
    * `cd experiments; aws s3 sync s3://thomas-moody/fmt/experiments/ .; cd ..`
    * `cd seq; aws s3 sync s3://thomas-moody/fmt/seq/ .; cd ..`

#### Building SINTAX DB
* working out of `seq` -> `cd seq`
* `mkdir db`
* Downloading species level RDP training set 18, formatted for UPARSE from [Edgar](https://drive5.com/usearch/manual/sintax_downloads.html)
* `cd db; wget https://drive5.com/sintax/rdp_16s_v18.fa.gz; cd ..;`
* extract gz :  `gunzip -c db/rdp_16s_v18.fa.gz > db/rdp_16s_v18.fasta`
* Build database : `usearch10 -makeudb_sintax db/rdp_16s_v18.fasta -output db/rdp_16s_v18.udb`

#### Continuing to Run 16S pipeline
* Should work w/ inc. memory
* need to pull down processing scripts
```bash
mkdir scripts
aws s3 cp s3://thomas-moody/fmt/2021-03-18_JxE_replication/scripts/process_otu.py scripts
aws s3 cp s3://thomas-moody/fmt/2021-03-18_JxE_replication/scripts/process_tax.py scripts
``` 
 
#### Complete run
* Succees! Important files
    * `cleaned_zotus.fa` -> Zotu sequences
    * `tax.csv` -> taxonomy table

### 2021-06-14
#### Rerunning w/ Redo Samples
* Needed to download and reorganize FMT03 reads, see that experiment folder for notes
* Cecum samples failed to sequence on first run
* They were re-seq'd on a second run
* Need to clean up dir first
* Remove pooled, derep, and compiled reads
    * `rm 01_pooled.fasta `
    * `rm 02_derep.fasta`
    * `rm 03_compiled.fastq`
* Remove tax table, zotu table, and zotu fasta's etc 
    * `rm tax.csv zmap.txt zotu.csv zotus.fa zotutab.txt`
    * `rm reads.sintax cleaned_zotus.fa`
