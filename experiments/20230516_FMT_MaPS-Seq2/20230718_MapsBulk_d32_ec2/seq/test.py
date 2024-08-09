import glob
import os


# Experiment list
# Add new experiemnts to this list to re-run pipeline
experiment_list = [
    "../experiments/2017-10-10_pairwise_cohousing",
    "../experiments/2018-09-30_JxE_FMT",
    "../experiments/2018-11-13_JxE_FMT02",
    "../experiments/2018-11-14_JxE_BCT01",
    "../experiments/2019-07-11_JxE_BCT02",
    "../experiments/2021-03-18_JxE_FMT03",
    "../experiments/2021-05-12_envigo_culture_pilot",
    "../experiments/2021-09-25_Ope_JvE_polysaccharide_screen"
]

# Collect seq_id's from filepaths
def collect_seq_ids():
    experiments = []
    seq_ids = []

    for experiment in experiment_list:
        reads = glob.glob("{}/seq/00_raw_reads/*.fastq.gz".format(experiment))
        for read in reads:
            seq_id = os.path.basename(read).split("_")[0]
            if seq_id not in seq_ids:
                seq_ids.append(seq_id)
                experiments.append(experiment)
    
    return seq_ids, experiments

seq_ids, experiments = collect_seq_ids()
print(seq_ids)