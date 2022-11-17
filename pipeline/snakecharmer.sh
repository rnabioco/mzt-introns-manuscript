#!/usr/bin/env bash
#BSUB -J mzintrons
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4] " 
#BSUB -q rna

mkdir -p logs

set -o nounset -o pipefail -o errexit -x

args=' -q rna 
       -o {log}.out 
       -e {log}.err 
       -J {rule} 
       -R "select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}] span[hosts=1] " 
       -n {threads}  ' 
    

# load modules to load singularity
. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load singularity

#### execute snakemake ####

run_snakemake () {
  snakemake \
      --drmaa "$args" \
      --snakefile Snakefile \
      --jobs 200 \
      --resources all_threads=200 \
      --latency-wait 50 \
      --printshellcmds \
      --rerun-incomplete  \
      --configfile $1 \
      --singularity-args '--bind /beevol/home/riemondy' \
      --use-singularity 
}

run_snakemake config-drosophila-rissland-v2.yaml 
run_snakemake config-drosophila-eichhorn-v2.yaml
run_snakemake config-xenopus-polya-v2.yaml
run_snakemake config-xenopus-ribozero-v2.yaml
run_snakemake config-zebrafish-v2.yaml 
run_snakemake config-coral.yaml
