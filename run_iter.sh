#!/bin/sh
source /broad/software/scripts/useuse
reuse -q R-4.0
cd /broad/thechenlab/Dylan/slideseq/spacexr/
set -e
mkdir -p logs
Rscript CANCER_de_script.R > logs/iter.txt