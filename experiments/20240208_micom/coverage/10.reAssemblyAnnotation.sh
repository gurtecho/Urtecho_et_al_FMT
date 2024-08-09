#!/usr/bin/env bash

chmod +x ./bin/*
mkdir -p output_09_reAssemblyAnnotation
for eachS in `cat ./supplier.list`
do
    echo $eachS
    checkm lineage_wf -t 96 --pplacer_threads 64 --tab_table -x fasta \
        -f ./output_09_reAssemblyAnnotation/$eachS\.checkM.stat \
        ./output_08_contigReAssemble/$eachS\_fasta \
        ./output_09_reAssemblyAnnotation/$eachS\.checkM

### GTdb-Tk cannot be ran in bash script and should be ran as command line
    gtdbtk classify_wf --genome_dir ./output_08_contigReAssemble/$eachS\_fasta \
        --out_dir ./output_09_reAssemblyAnnotation/$eachS\.gtdbtk \
        -x fasta --prefix $eachS --cpus 45
    cp ./output_09_reAssemblyAnnotation/$eachS\.gtdbtk/$eachS\.bac120.summary.tsv \
        ./output_09_reAssemblyAnnotation/$eachS\.gtdbtk.stat 

    tar czvf ./output_09_reAssemblyAnnotation/$eachS\.checkM.tar.gz ./output_09_reAssemblyAnnotation/$eachS\.checkM &
    tar czvf ././output_09_reAssemblyAnnotation/$eachS\.gtdbtk.tar.gz ././output_09_reAssemblyAnnotation/$eachS\.gtdbtk &

    ls -l output_08_contigReAssemble/$eachS\_fasta  | awk '{print "./output_08_contigReAssemble/'$eachS'_fasta/"$9}' > output_08_contigReAssemble/$eachS\_fasta.list

    for eachS2 in `cat ./supplier.list`
    do
        fastANI --rl ./output_08_contigReAssemble/$eachS\_fasta.list \
                --ql ./output_08_contigReAssemble/$eachS2\_fasta.list \
                -t 32 -o ./output_09_reAssemblyAnnotation/fastANI.$eachS\.$eachS2\.stat
    done
done

wait


