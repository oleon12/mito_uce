ls -lh results/S_bogotensis_AMNH_207854/bam/

CONDA_ENV="mt_pipeline"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

for bam in results/*/bam/*.mapped.bam; do
    sample=$(basename "$bam" .mapped.bam)
    mapped=$(samtools view -c "$bam")
    echo -e "$sample\t$mapped"
done
