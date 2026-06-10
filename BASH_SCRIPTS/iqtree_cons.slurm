#!/bin/bash
#SBATCH --job-name=iqtree_cons
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/iqtree_cons_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/iqtree_cons_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/odl7/sturnira_mito

source ~/miniforge3/etc/profile.d/conda.sh
conda activate iqtree

ALIGN="results/alignments/combined_masked_consensus_intersection.mafft.fasta"
OUTDIR="results/iqtree_consensus"
OUTGROUP="Artibeus_PP853570.1,Glossophaga_NC_065682.1"

mkdir -p "$OUTDIR"

if [[ ! -f "$ALIGN" ]]; then
    echo "ERROR: alignment not found: $ALIGN"
    exit 1
fi

iqtree \
  -s "$ALIGN" \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T AUTO \
  -o "$OUTGROUP" \
  --prefix "$OUTDIR/combined_masked_consensus"

echo "Done. IQ-TREE consensus results in: $OUTDIR"
