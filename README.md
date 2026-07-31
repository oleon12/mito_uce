# Mito-UCE

Mito-UCE is a reproducible workflow for recovering mitochondrial sequence data from paired-end UCE reads by mapping reads to a complete mitochondrial reference. The workflow produces two complementary datasets:

1. A **primary 13-protein-coding-gene (13-PCG) dataset**, aligned gene by gene through amino-acid translation and codon-aware back-translation.
2. An **optional complete-mitogenome dataset**, retained as a complementary sensitivity analysis after alignment and manual quality review.

This guide is based on and expands the original [mtDNA-MitoFinder pipeline](https://github.com/Agustol/mtDNA-mitofinder-pipeline). Credit for the original workflow goes to [Agusto Luzuriaga-Neira](https://github.com/Agustol).

> [!IMPORTANT]
> The current repository prepares, validates, filters, concatenates, and aligns mitochondrial sequences. It does **not** prescribe or include a final IQ-TREE analysis. Use the resulting alignments and partition files with the phylogenetic software and inference strategy appropriate for your study.

---

## Recommended directory structure

Organizing the project as shown below allows the scripts to run with minimal path changes.

```text
main/
├── BASH_SCRIPTS/
├── CONFS/
│   ├── environment.yml
│   ├── mitofinder_env.yml
│   ├── sample_list.txt
│   ├── out_group_template.txt
│   └── outgroup_list.txt          # created from the template
├── ERR_OUT/
├── LOG/
├── raw_data/
│   ├── outgroups/
│   └── UCE/
├── references/
└── results/
```

Most scripts currently use:

```bash
cd /scratch/odl7/sturnira_mito
```

Change `WORKDIR` or the initial `cd` command in every script when adapting the workflow to another project.

---

## Workflow overview

```text
Reference FASTA + GenBank
        │
        ├── prep_ref.slurm
        └── cds_bed.slurm
                 │
Paired UCE FASTQs
        │
        └── map_all_ref.slurm
                 │
        duplicate-marked mapped BAMs
                 │
        └── vcf_file.slurm
                 │
        normalized haploid audit VCFs
                 │
        ├── make_consensus.slurm            optional diagnostic output
        └── make_masked_consensus.slurm     required publication output
                 │
        ├── complete masked mitogenomes
        └── extract_cds.slurm
                 │
        13 named PCGs per specimen
                 │
Public outgroup FASTA + GenBank
        │
        └── outgroup_from_gb.sh
                 │
Independent M40 / M30 / M20 filtering
        │
        ├── concat_cds.sh
        └── concat_cons.sh
                 │
        ├── mafft_cds.slurm
        └── mafft_cons.slurm
                 │
Validated alignments for downstream phylogenetic analysis
```

---

# 0. Software environments

The main mapping and consensus workflow uses the Conda environment defined in [`CONFS/environment.yml`](https://github.com/oleon12/mito_uce/blob/main/CONFS/environment.yml).

```bash
conda env create -f CONFS/environment.yml
```

The scripts expect this environment to be named:

```bash
mt_pipeline
```

The repository also retains [`CONFS/mitofinder_env.yml`](https://github.com/oleon12/mito_uce/blob/main/CONFS/mitofinder_env.yml) for MITOFinder-compatible or legacy analyses:

```bash
conda env create -f CONFS/mitofinder_env.yml --solver=classic
```

The current reference-mapping workflow documented below uses `mt_pipeline`.

Core software includes:

- BWA
- SAMtools
- BCFtools
- MAFFT
- Biopython
- standard Unix tools

---

# 1. Input configuration

## 1.1. Ingroup sample list

Create [`CONFS/sample_list.txt`](https://github.com/oleon12/mito_uce/blob/main/CONFS/sample_list.txt) with exactly three whitespace-separated fields per non-comment line:

```text
SAMPLE_ID R1_FASTQ R2_FASTQ
```

Example:

```text
S_angeli_AMNH_214197 raw_data/UCE/CBT03_L0074_R1.fastq.gz raw_data/UCE/CBT03_L0074_R2.fastq.gz
S_bogotensis_AMNH_246573 raw_data/UCE/CBT03_L0096_R1.fastq.gz raw_data/UCE/CBT03_L0096_R2.fastq.gz
```

Sample identifiers should contain only letters, numbers, periods, underscores, and hyphens. Each identifier must be unique because it becomes the BAM read-group sample name, VCF sample name, FASTA identifier, and downstream taxon name.

## 1.2. Outgroup manifest

The repository includes [`CONFS/out_group_template.txt`](https://github.com/oleon12/mito_uce/blob/main/CONFS/out_group_template.txt). Copy it to the filename expected by the outgroup-preparation script:

```bash
cp CONFS/out_group_template.txt CONFS/outgroup_list.txt
```

Then edit `CONFS/outgroup_list.txt`. Each non-comment line must contain:

```text
OUTGROUP_ID FASTA_PATH GENBANK_PATH
```

Example:

```text
Artibeus_PP853570.1 raw_data/outgroups/Artibeus_PP853570.1.fasta raw_data/outgroups/Artibeus_PP853570.1.gb
Glossophaga_NC_065682.1 raw_data/outgroups/Glossophaga_NC_065682.1.fasta raw_data/outgroups/Glossophaga_NC_065682.1.gb
```

The FASTA and GenBank files for each outgroup must represent the same complete circular mitochondrial molecule.

---

# 2. Reference genome

A complete mitochondrial reference is required in both FASTA and GenBank formats. Both files must contain exactly one record and must represent the same sequence.

The example project uses:

```text
references/S_ludovici_QCAZ_18312.fasta
references/S_ludovici_QCAZ_18312.gb
```

## 2.1. Prepare and validate the reference

Run [`prep_ref.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/prep_ref.slurm).

```bash
sbatch BASH_SCRIPTS/prep_ref.slurm
```

The script:

- validates exact FASTA–GenBank sequence identity;
- requires one unambiguous circular DNA sequence;
- validates the conventional mammalian tRNA-Phe circular origin;
- verifies all 13 mitochondrial protein-coding genes;
- checks `codon_start`, vertebrate mitochondrial translation table 2, annotated translations, and internal stops;
- builds and validates all classic BWA index components;
- creates and validates the SAMtools `.fai` index;
- detects stale indexes through checksums;
- performs functional BWA-MEM and `samtools faidx` smoke tests.

Essential outputs include:

```text
references/S_ludovici_QCAZ_18312.fasta.amb
references/S_ludovici_QCAZ_18312.fasta.ann
references/S_ludovici_QCAZ_18312.fasta.bwt
references/S_ludovici_QCAZ_18312.fasta.pac
references/S_ludovici_QCAZ_18312.fasta.sa
references/S_ludovici_QCAZ_18312.fasta.fai

results/reference_preparation/reference_qc.tsv
results/reference_preparation/reference_cds_qc.tsv
results/reference_preparation/reference_index_manifest.tsv
```

## 2.2. Build the reference CDS coordinates

Run [`cds_bed.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/cds_bed.slurm).

```bash
sbatch BASH_SCRIPTS/cds_bed.slurm
```

The script parses the GenBank annotation structurally with Biopython. It validates and extracts the coding portions of:

```text
ND1 ND2 COX1 COX2 ATP8 ATP6 COX3 ND3 ND4L ND4 ND5 ND6 CYTB
```

It preserves strand and compound-feature information, excludes annotated terminal stop bases, validates translations, and creates a coding-only BED12 file.

Essential outputs:

```text
results/cds_coords.bed
results/cds_coords.tsv
```

Do not replace this BED12 file with a simple six-column BED. The revised extraction workflow uses its block and strand information.

---

# 3. Map UCE reads to the mitochondrial reference

Run [`map_all_ref.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/map_all_ref.slurm).

```bash
sbatch BASH_SCRIPTS/map_all_ref.slurm
```

For every entry in `sample_list.txt`, the script:

- validates paired FASTQ structure and mate names;
- adds an `@RG` read group with `SM=SAMPLE_ID`;
- maps with BWA-MEM;
- retains primary, mapped, QC-passed alignments;
- name-sorts and runs `samtools fixmate`;
- coordinate-sorts the BAM;
- marks, but does not delete, duplicate reads;
- validates BAM integrity, sample identity, flags, and indexing;
- finalizes outputs only after all checks pass.

The primary BAM output is:

```text
results/SAMPLE/bam/SAMPLE.mapped.bam
results/SAMPLE/bam/SAMPLE.mapped.bam.bai
```

The workflow no longer requires a large intermediate `SAMPLE.sorted.bam`.

Additional audit files include:

```text
results/SAMPLE/bam/SAMPLE.mapping_qc.tsv
results/SAMPLE/bam/SAMPLE.bwa_mem.log
results/SAMPLE/bam/SAMPLE.markdup_stats.txt
results/SAMPLE/bam/SAMPLE.mapping_commands.log
```

## 3.1. Mapping summary

Run [`mapping_summary.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/mapping_summary.slurm).

```bash
sbatch BASH_SCRIPTS/mapping_summary.slurm
```

The summary uses the same quality-filtered evidence later used for masking:

- minimum base quality 20;
- minimum mapping quality 20;
- duplicate reads excluded;
- overlapping paired-read segments counted once;
- zero-depth positions included.

Output:

```text
results/mapping_summary.tsv
results/mapping_summary_failed_samples.txt
```

Reported metrics include mapped and duplicate counts, mean and median genome-wide depth, maximum depth, coverage at ≥1×, ≥3×, ≥5×, and ≥10×, and a heuristic warning for reduced coverage near the artificial boundary of the linearized circular genome.

## 3.2. Test mapping outputs

Run [`test_map_ref.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/test_map_ref.sh).

```bash
bash BASH_SCRIPTS/test_map_ref.sh
```

This is a read-only regression test. It validates BAM integrity, indexes, read groups, duplicate marking, alignment flags, reference identity, and agreement with `mapping_summary.tsv`.

Output:

```text
results/tests/test_map_ref.tsv
```

---

# 4. Call mitochondrial variants

Run [`vcf_file.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/vcf_file.slurm).

```bash
sbatch BASH_SCRIPTS/vcf_file.slurm
```

The script creates one normalized, haploid, soft-filtered audit VCF per sample. It:

- requires the revised duplicate-marked BAM and matching read-group sample name;
- uses explicit mapping-quality and base-quality thresholds;
- chooses a sample-specific mpileup depth cap above the observed maximum depth;
- applies full BAQ;
- calls haploid genotypes;
- normalizes records against the reference;
- splits multiallelic records;
- retains SNPs, indels, complex calls, and low-confidence calls with informative FILTER labels.

The final VCF is:

```text
results/SAMPLE/vcf/SAMPLE.vcf.gz
results/SAMPLE/vcf/SAMPLE.vcf.gz.csi
```

A record is eligible for later consensus application only when it is a high-confidence haploid ALT SNP that passes the configured QUAL, depth, alternate-depth, alternate-fraction, and proximity-to-indel rules.

Other records remain in the VCF so that the masked-consensus step can treat their reference spans as uncertain rather than silently retaining the reference allele.

Additional outputs:

```text
results/SAMPLE/vcf/SAMPLE.raw.normalized.vcf.gz
results/SAMPLE/vcf/SAMPLE.variant_qc.tsv
results/SAMPLE/vcf/SAMPLE.bcftools_stats.txt
results/SAMPLE/vcf/SAMPLE.vcf_commands.log
results/vcf_summary.tsv
results/vcf_failed_samples.txt
```

> [!NOTE]
> This workflow builds a conservative homoplasmic consensus. It is not a dedicated heteroplasmy analysis.

## 4.1. Test VCF outputs

Run [`test_vcf.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/test_vcf.sh).

```bash
bash BASH_SCRIPTS/test_vcf.sh
```

The test validates:

- VCF and CSI integrity;
- sample identity;
- haploid genotypes;
- reference coordinates;
- FILTER definitions;
- PASS SNP criteria;
- agreement with `vcf_summary.tsv`.

Output:

```text
results/tests/test_vcf.tsv
```

---

# 5. Build mitochondrial consensuses

## 5.1. Optional unmasked diagnostic consensus

[`make_consensus.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/make_consensus.slurm) produces an unmasked, SNP-only diagnostic reconstruction.

```bash
sbatch BASH_SCRIPTS/make_consensus.slurm
```

Output:

```text
results/SAMPLE/consensus_diagnostic/SAMPLE_mt_unmasked_snp_only.fasta
results/unmasked_consensus_summary.tsv
```

> [!WARNING]
> These sequences retain reference bases at positions with insufficient or absent evidence. They are for diagnostics only and must not be used as phylogenetic input.

The optional diagnostic regression test is:

```bash
bash BASH_SCRIPTS/test_consensus.sh
```

Output:

```text
results/tests/test_consensus.tsv
```

## 5.2. Required coordinate-preserving masked consensus

Run [`make_masked_consensus.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/make_masked_consensus.slurm).

```bash
sbatch BASH_SCRIPTS/make_masked_consensus.slurm
```

This is the required consensus step for all downstream analyses.

The script:

- applies only high-confidence PASS haploid SNPs;
- never inserts or deletes sequence relative to the reference;
- masks called indels, complex records, low-confidence calls, and uncertain genotypes;
- masks positions below the configured quality-filtered depth;
- calculates depth after base-quality and mapping-quality filtering;
- counts overlapping paired-read segments once;
- requires the final consensus to equal the reference length exactly;
- writes the correct `>SAMPLE` FASTA header;
- builds and validates the FASTA index.

The default depth threshold is:

```bash
MIN_DEPTH=3
```

A separate 5× sensitivity reconstruction may be useful for especially incomplete or unstable specimens.

Essential outputs:

```text
results/SAMPLE/consensus_masked/SAMPLE_mt_masked.fasta
results/SAMPLE/consensus_masked/SAMPLE_mt_masked.fasta.fai
results/SAMPLE/consensus_masked/SAMPLE.masked_consensus_qc.tsv

results/masked_consensus_summary.tsv
results/masked_consensus_failed_samples.txt
```

Per-sample audit outputs also preserve depth tables, mask BED files, normalized/accepted VCF subsets, and command logs.

## 5.3. Test masked consensuses

Run [`test_masked_consensus.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/test_masked_consensus.sh).

```bash
bash BASH_SCRIPTS/test_masked_consensus.sh
```

The test validates sequence identifiers, symbols, exact reference length, `.fai` consistency, full-sequence retrieval, missing-data counts, and agreement with `masked_consensus_summary.tsv`.

Output:

```text
results/tests/test_masked_consensus.tsv
```

---

# 6. Extract the 13 mitochondrial PCGs

Run [`extract_cds.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/extract_cds.slurm).

```bash
sbatch BASH_SCRIPTS/extract_cds.slurm
```

The script extracts the 13 coding regions from each fixed-coordinate masked consensus using `results/cds_coords.bed` and `results/cds_coords.tsv`.

Each sample produces exactly 13 records in biological gene order:

```fasta
>SAMPLE|ND1
>SAMPLE|ND2
>SAMPLE|COX1
...
>SAMPLE|CYTB
```

Output:

```text
results/SAMPLE/genes_masked/SAMPLE_13PCG_masked.fasta
results/SAMPLE/genes_masked/SAMPLE_13PCG_qc.tsv

results/extract_masked_cds_summary.tsv
results/extract_masked_cds_failed_samples.txt
```

The script validates:

- the exact 13-gene set and order;
- valid IUPAC nucleotides;
- expected ingroup gene lengths;
- coding frame;
- mitochondrial translation table;
- absence of definite internal stop codons.

## 6.1. Test extracted CDS data

Run [`test_extract_cds.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/test_extract_cds.sh).

```bash
bash BASH_SCRIPTS/test_extract_cds.sh
```

The test independently validates ingroup and outgroup 13-PCG FASTAs, gene names and order, frames, stop codons, coding lengths, missing-data counts, and agreement with the production summaries.

Output:

```text
results/tests/test_extract_cds.tsv
```

---

# 7. Prepare outgroups

Place each outgroup FASTA and matching GenBank file under:

```text
raw_data/outgroups/
```

Create `CONFS/outgroup_list.txt` from [`CONFS/out_group_template.txt`](https://github.com/oleon12/mito_uce/blob/main/CONFS/out_group_template.txt), as described in Section 1.2.

Then run [`outgroup_from_gb.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/outgroup_from_gb.sh).

```bash
bash BASH_SCRIPTS/outgroup_from_gb.sh
```

For every outgroup, the script:

- verifies that FASTA and GenBank represent the same circular sequence;
- validates all 13 PCG annotations and translations;
- extracts every gene from that outgroup's own GenBank annotation;
- writes named, frame-valid coding sequences;
- standardizes complete-genome orientation and circular origin to the reference;
- validates PCG order and strand pattern;
- writes complete QC and audit reports.

Outputs include:

```text
results/OUTGROUP/consensus_masked/OUTGROUP_mt_masked.fasta
results/OUTGROUP/consensus_masked/OUTGROUP_mt_masked.fasta.fai

results/OUTGROUP/genes_masked/OUTGROUP_13PCG_masked.fasta
results/OUTGROUP/genes_masked/OUTGROUP_13PCG_qc.tsv

results/OUTGROUP/outgroup_qc/OUTGROUP_circular_standardization.tsv

results/outgroup_samples.txt
results/outgroup_extraction_summary.tsv
results/outgroup_extraction_failed.txt
```

The `consensus_masked` directory name is retained for compatibility with the ingroup workflow. Public outgroup genomes are complete standardized sequences, not reference-derived masked consensuses.

---

# 8. Filter samples by missing data

The workflow creates three nested ingroup datasets:

```text
M40: percent N ≤40%
M30: percent N ≤30%
M20: percent N ≤20%
```

These thresholds allow explicit sensitivity analyses to missing data.

## 8.1. Filter the 13-PCG data

Run [`filter_from_cds.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/filter_from_cds.sh).

```bash
bash BASH_SCRIPTS/filter_from_cds.sh
```

The script reads:

```text
results/extract_masked_cds_summary.tsv
results/cds_coords.tsv
CONFS/sample_list.txt
```

It verifies summary structure and recalculates missing-data percentages before writing:

```text
results/keep_samples_cds_le40.txt
results/keep_samples_cds_le30.txt
results/keep_samples_cds_le20.txt

results/drop_samples_cds_gt40.txt
results/drop_samples_cds_gt30.txt
results/drop_samples_cds_gt20.txt

results/filter_from_cds_report.tsv
results/filter_from_cds_counts.tsv
```

## 8.2. Filter complete masked mitogenomes

Run [`filter_from_con.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/filter_from_con.sh).

```bash
bash BASH_SCRIPTS/filter_from_con.sh
```

The script reads and verifies the actual FASTAs against:

```text
results/masked_consensus_summary.tsv
```

It writes:

```text
results/keep_samples_cons_le40.txt
results/keep_samples_cons_le30.txt
results/keep_samples_cons_le20.txt

results/drop_samples_cons_gt40.txt
results/drop_samples_cons_gt30.txt
results/drop_samples_cons_gt20.txt

results/filter_from_con_report.tsv
results/filter_from_con_counts.tsv
```

## 8.3. Optional matched-taxon intersections

Run [`filter_intersection.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/filter_intersection.sh) only when you need a matched-taxon comparison between the 13-PCG and complete-mitogenome analyses.

```bash
bash BASH_SCRIPTS/filter_intersection.sh
```

Outputs:

```text
results/keep_samples_intersection_le40.txt
results/keep_samples_intersection_le30.txt
results/keep_samples_intersection_le20.txt
results/filter_intersection_report.tsv
results/filter_intersection_counts.tsv
```

> [!IMPORTANT]
> Do not use the intersection lists as the default filter for the primary analyses. The 13-PCG and complete-mitogenome datasets should normally use their independent CDS and consensus keep lists. The intersection is an optional sensitivity analysis that holds taxon sampling constant.

---

# 9. Build unaligned analysis matrices

## 9.1. Prepare gene-wise 13-PCG matrices

Run [`concat_cds.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/concat_cds.sh).

```bash
bash BASH_SCRIPTS/concat_cds.sh
```

Despite its historical filename, this script does **not** immediately concatenate unaligned genes into a supermatrix. It validates every taxon and creates one unaligned multi-FASTA per gene for M40, M30, and M20.

Outputs:

```text
results/cds_gene_matrices/
├── M40/
│   ├── taxa.txt
│   └── unaligned_by_gene/
│       ├── ND1.fasta
│       ├── ND2.fasta
│       ├── COX1.fasta
│       └── ...
├── M30/
└── M20/
```

The script includes the configured outgroups and validates complete taxon and gene sets before finalizing the matrices.

## 9.2. Prepare complete-mitogenome matrices

Run [`concat_cons.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/concat_cons.sh).

```bash
bash BASH_SCRIPTS/concat_cons.sh
```

This creates unaligned complete-mitogenome multi-FASTA datasets for the independently filtered M40, M30, and M20 consensus sets.

Outputs:

```text
results/consensus_matrices/M40/Sturnira_complete_mt_M40.unaligned.fasta
results/consensus_matrices/M30/Sturnira_complete_mt_M30.unaligned.fasta
results/consensus_matrices/M20/Sturnira_complete_mt_M20.unaligned.fasta
```

The complete-mitogenome analysis is complementary to, not independent from, the 13-PCG analysis because it contains the same protein-coding genes plus rRNAs, tRNAs, and noncoding regions.

---

# 10. Align the matrices

## 10.1. Protein-guided, codon-aware 13-PCG alignment

Run [`mafft_cds.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/mafft_cds.slurm).

```bash
sbatch BASH_SCRIPTS/mafft_cds.slurm
```

This SLURM array processes M40, M30, and M20. For every gene, it:

1. translates nucleotide CDS under vertebrate mitochondrial code 2;
2. aligns amino-acid sequences with MAFFT;
3. back-translates the protein alignment into codons;
4. verifies codon gaps, frames, lengths, taxa, and stop codons;
5. concatenates the 13 aligned genes in fixed biological order;
6. writes a 39-partition gene-by-codon-position file.

Example M40 outputs:

```text
results/alignments/cds/M40/
├── protein_inputs/
├── protein_alignments/
├── codon_alignments/
├── mafft_logs/
├── Sturnira_13PCG_M40.codon_aligned.fasta
├── Sturnira_13PCG_M40.partitions.txt
├── Sturnira_13PCG_M40.coordinates.tsv
├── Sturnira_13PCG_M40.alignment_qc.tsv
├── coding_overlap_report.tsv
├── alignment_strategy.tsv
└── mafft_commands.tsv
```

The overlap report documents mitochondrial coding overlaps that are represented in more than one complete gene.

## 10.2. Complete-mitogenome alignment

Run [`mafft_cons.slurm`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/mafft_cons.slurm).

```bash
sbatch BASH_SCRIPTS/mafft_cons.slurm
```

This SLURM array aligns the M40, M30, and M20 complete-mitogenome matrices. It:

- temporarily adds the reference as an alignment anchor;
- retains an alignment containing the reference for annotation mapping;
- removes the temporary anchor from the final phylogenetic matrix;
- verifies that each ungapped sequence is unchanged;
- writes a reference-coordinate-to-alignment-column map;
- reports per-taxon and per-column missingness and gaps;
- flags possible circular-origin or edge problems.

Example M40 outputs:

```text
results/alignments/consensus/M40/
├── Sturnira_complete_mt_M40.mafft.fasta
├── Sturnira_complete_mt_M40.mafft.with_reference.fasta
├── Sturnira_complete_mt_M40.reference_coordinate_map.tsv
├── Sturnira_complete_mt_M40.alignment_qc.tsv
├── Sturnira_complete_mt_M40.column_stats.tsv
├── alignment_strategy.tsv
└── mafft.log
```

No automatic trimming is performed. The complete-mitogenome alignment, especially the control region and other difficult noncoding regions, must be reviewed before phylogenetic inference.

---

# 11. Downstream phylogenetic inference

The repository intentionally stops after validated alignment generation.

For the primary 13-PCG analysis, use:

```text
results/alignments/cds/M40/Sturnira_13PCG_M40.codon_aligned.fasta
results/alignments/cds/M40/Sturnira_13PCG_M40.partitions.txt
```

with corresponding M30 and M20 files for missing-data sensitivity analyses.

For the optional complete-mitogenome analysis, use the reviewed matrices under:

```text
results/alignments/consensus/M40/
results/alignments/consensus/M30/
results/alignments/consensus/M20/
```

Choose and document the phylogenetic software, substitution-model strategy, partition treatment, support metrics, random seeds, and convergence or stability checks appropriate for your research question.

---

# 12. Retired scripts

The following historical scripts are retired:

- [`rename_fasta.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/rename_fasta.sh)
- [`rename_fai.sh`](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/rename_fai.sh)

Do not run them.

All current production scripts write correct FASTA identifiers and indexes when the files are created. Post hoc renaming or routine reindexing could invalidate QC reports, conceal altered files, or destroy biological gene names.

When a validated FASTA or index is missing or inconsistent, rerun the script that produced it.

---

# 13. Recommended execution order

```bash
# 1. Reference
sbatch BASH_SCRIPTS/prep_ref.slurm
sbatch BASH_SCRIPTS/cds_bed.slurm

# 2. Mapping
sbatch BASH_SCRIPTS/map_all_ref.slurm
sbatch BASH_SCRIPTS/mapping_summary.slurm
bash BASH_SCRIPTS/test_map_ref.sh

# 3. Variant calling
sbatch BASH_SCRIPTS/vcf_file.slurm
bash BASH_SCRIPTS/test_vcf.sh

# 4. Required masked consensus
sbatch BASH_SCRIPTS/make_masked_consensus.slurm
bash BASH_SCRIPTS/test_masked_consensus.sh

# 5. Extract coding genes
sbatch BASH_SCRIPTS/extract_cds.slurm

# 6. Prepare public outgroups
bash BASH_SCRIPTS/outgroup_from_gb.sh

# 7. Validate all 13-PCG FASTAs
bash BASH_SCRIPTS/test_extract_cds.sh

# 8. Independent missing-data filters
bash BASH_SCRIPTS/filter_from_cds.sh
bash BASH_SCRIPTS/filter_from_con.sh

# Optional matched-taxon sensitivity lists
bash BASH_SCRIPTS/filter_intersection.sh

# 9. Build matrices
bash BASH_SCRIPTS/concat_cds.sh
bash BASH_SCRIPTS/concat_cons.sh

# 10. Align M40, M30, and M20
sbatch BASH_SCRIPTS/mafft_cds.slurm
sbatch BASH_SCRIPTS/mafft_cons.slurm
```

Optional diagnostic consensus:

```bash
sbatch BASH_SCRIPTS/make_consensus.slurm
bash BASH_SCRIPTS/test_consensus.sh
```

---

# 14. Primary reports to inspect

Before downstream phylogenetic inference, inspect at least:

```text
results/reference_preparation/reference_qc.tsv
results/mapping_summary.tsv
results/vcf_summary.tsv
results/masked_consensus_summary.tsv
results/extract_masked_cds_summary.tsv
results/outgroup_extraction_summary.tsv
results/filter_from_cds_report.tsv
results/filter_from_con_report.tsv
results/tests/test_map_ref.tsv
results/tests/test_vcf.tsv
results/tests/test_masked_consensus.tsv
results/tests/test_extract_cds.tsv
```

Every required production and regression-test row should have:

```text
status = passed
```

The M40, M30, and M20 trees or other downstream results should then be compared explicitly to evaluate sensitivity to missing data and taxon inclusion.
