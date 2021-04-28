#!/bin/zsh

source ~/.zprofile
# bash wrapper script to plot karyogram
conda activate py36-plot
python ~/sherlock/oak/LocalAncestry/scripts/plotting/plot_karyogram.py \
    --bed_path ~/sherlock/oak/LocalAncestry/data/bed/chr1/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr1/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr2/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr2/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr3/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr3/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr4/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr4/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr5/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr5/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr6/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr6/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr7/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr7/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr8/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr8/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr9/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr9/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr10/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr10/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr11/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr11/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr12/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr12/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr13/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr13/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr14/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr14/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr15/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr15/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr16/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr16/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr17/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr17/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr18/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr18/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr19/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr19/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr20/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr20/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr21/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr21/NWD149502.B.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr22/NWD149502.A.bed \
               ~/sherlock/oak/LocalAncestry/data/bed/chr22/NWD149502.B.bed \
    --ind NWD149502 \
    --centromeres ~/sherlock/oak/LocalAncestry/data/maps/centromere.tsv \
    --pop_order CEU,YRI \
    --colors 455A8E,E48671,22223A \
    --out ~/sherlock/oak/LocalAncestry/results/v2/plots/NWD149502.png
conda deactivate
