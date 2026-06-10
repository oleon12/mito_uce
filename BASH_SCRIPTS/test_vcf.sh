CONDA_ENV="mt_pipeline"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"


for vcf in results/*/vcf/*.vcf.gz; do
    sample=$(basename "$vcf" .vcf.gz)
    nvar=$(bcftools view -H "$vcf" | wc -l)
    echo -e "$sample\t$nvar"
done
