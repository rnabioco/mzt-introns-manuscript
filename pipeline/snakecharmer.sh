#!/usr/bin/env bash
#BSUB -J mztintrons
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4] " 
#BSUB -q rna

mkdir -p logs

set -o nounset -o pipefail -o errexit -x

# cluster specific arguments
args=' -q rna 
       -o {log}.out 
       -e {log}.err 
       -J {rule} 
       -R "select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}] span[hosts=1] " 
       -n {threads}  ' 

#### execute snakemake ####

snakemake --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 15 \
    --resources all_threads=7 \
    --latency-wait 50 \
    --printshellcmds \
    --rerun-incomplete  \
    --configfile config-test.yaml 
