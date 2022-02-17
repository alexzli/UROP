#!/bin/sh
source /broad/software/scripts/useuse
reuse -q R-4.0
cd /broad/thechenlab/Alex/UROP/scripts
set -e
mkdir -p logs
Rscript gen_simulated_script.R > logs/iter.txt