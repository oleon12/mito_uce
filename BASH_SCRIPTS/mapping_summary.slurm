#!/bin/bash
#SBATCH --job-name=summ_mapping
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=244G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/summ_mapping_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/summ_mapping_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail


######

cd /scratch/odl7/sturnira_mito

######

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mt_pipeline


OUTDIR="results"
SUMMARY="results/mapping_summary.tsv"

OUTDIR=$(realpath -m "$OUTDIR")
SUMMARY=$(realpath -m "$SUMMARY")

mkdir -p "$(dirname "$SUMMARY")"

echo -e "Sample\tMapped_Reads\tCovered_Positions\tAverage_Depth" > "$SUMMARY"

for bam in "$OUTDIR"/*/bam/*.mapped.bam; do
    [[ -f "$bam" ]] || continue

    sample=$(basename "$bam" .mapped.bam)
    mapped=$(samtools view -c "$bam")

    depth_stats=$(samtools depth "$bam" | awk '
        {sum+=$3; n++}
        END{
            if(n>0){
                printf "%d\t%.4f\n", n, sum/n
            } else {
                printf "0\t0\n"
            }
        }')

    covered=$(echo "$depth_stats" | cut -f1)
    avgdepth=$(echo "$depth_stats" | cut -f2)

    echo -e "${sample}\t${mapped}\t${covered}\t${avgdepth}" >> "$SUMMARY"
done

echo "Summary written to: $SUMMARY"
