
## Organizing Experiment Data
* Working on EC2 to pull down
* Working off YH S3 master metadata sheet [mice_16SV4_master.xlsx](s3://yiming-huang/20201217_16SV4_mouseFMT_master/mice_16SV4_master.xlsx)
* Generated a metadata sheet based on above : [samples.csv](metadata/samples.csv)
* Original Sequencing data on RS S3 : [20171010_16S_FMT_heatmap](s3://ravi_sheth/seq_files/20171010_16S_FMT_heatmap/)
* Reads are stored in YH S3 : [rawdata](s3://yiming-huang/20201217_16SV4_mouseFMT_master/rawdata/)

## Uploading to EC2
* move existing metadata : `scp -i ~/Dropbox/aws_pems/tum.pem -r ../2017-10-10_pairwise_cohousing ec2-user@52.72.35.224:~/fmt/experiments`