#!/bin/bash
#SBATCH --job-name=reindex_cons
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/reindex_cons_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/reindex_cons_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/odl7/sturnira_mito

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mt_pipeline

for fasta in results/*/consensus_masked/*_mt_masked.fasta; do
    [[ -f "$fasta" ]] || continue
    rm -f "${fasta}.fai"
    samtools faidx "$fasta"
    echo "Reindexed: $fasta"
done
