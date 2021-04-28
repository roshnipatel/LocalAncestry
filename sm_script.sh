#!/usr/bin/bash

# Prefer to run this in an interactive session so you can see the jobs running
# and when they fail - don't really see a reason to run it as sbatch

snakemake --rerun-incomplete --keep-going -j 20 \
        --latency-wait 60 --wait-for-files \
        --cluster-config scripts/sm_slurm_config.json \
        --cluster "sbatch -p {cluster.queue} \
                        -t {cluster.time} \
                        --ntasks-per-node={cluster.tasks} \
                        --job-name={cluster.name} \
                        -o {cluster.output} \
                        -e {cluster.error} \
                        --nodes={cluster.nodes} \
                        --mem={cluster.memory}"
