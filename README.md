# Mito - UCE

This is a brief guide or tutorial on obtaining mitochondrial DNA, whether in the form of complete genomes or genes from UCE data. For this guide to work, your UCE data must be enriched during sequencing. This guide is a modification of the original [pipeline](https://github.com/Agustol/mtDNA-mitofinder-pipeline), so all credits go to [Agusto Luzuriaga-Neira](https://github.com/Agustol)

Also, to reduce mistakes, I strongly suggest organizing your data in the following structure so you do not need to change many of the scripts.

```
main/
  ├── BASH_SCRIPTS
  ├── CONFS
  ├── ERR_OUT
  ├── LOG
  ├── raw_data/
  │   ├── outgroups/
  │   └── UCE/
  ├── references/
  └── results/
```

## 0. Environments

Before you start the guide, you must build specific environments. In fact, you must create two different environments, one for samtools, bwa, beftools, bedtools, and another for Python 2.7, blast, spades, and mitofinder-compatible packages. You need to use conda to build the two environment (.yml) files in the <b>CONFS</b> folder.

```
# Make sure that you are in the main before run these lines.

conda env create -f CONFS/environment.yml
conda env create -f CONFS/mitofinder_env.yml --solver=classic

```

## 1. Reference genome

To extract the genome or genes from your data, you need a reference genome. In this case, I was working with species of the genus Sturnira, so I went to GenBank and downloaded the .fasta and .gb files for a complete genome of <b><i>Sturnira ludovici</i></b>. Thus, you will need to find a reference genome and save it into the <b>references</b> folder.

### 1.1. Prepare reference 

The first step is to prepare the reference files from the reference genome. For this step, you will run the script <b>prep_ref.slurm</b>. This script first builds the BWA index, which creates several files: .bwt, .sa, .amb, .ann, and .pac. Then it will create a FASTA index with samtools, producing a .fai file (tab-separated: chromosome name, length, offset, etc.). For this script, you must make sure that you call the proper environment with conda, set the path to the reference fasta file, and set the proper working directory (main).

All files will be saved in the <b>references</b> folder.

```
# Set working directory
cd /scratch/odl7/sturnira_mito

# Call environment
CONDA_ENV="mt_pipeline"

# Set references path
REF="references/S_ludovici_QCAZ_18312.fasta"

```

### 1.2. Build CDS beds

Now, you must run the script <b>cds_bed.slurm</b>. This script takes a GenBank (.gb) file and a matching FASTA file, extracts all CDS (coding sequence) coordinates, and converts them into BED format (0‑based, half‑open) for downstream analyses. For this script, you make sure to set the working directory, the fasta file, and the gb file of the reference genome. The output files will be saved in the <b>results</b> folder. <i>You can change the output file name if you want, but you will need to update it in subsequent steps.</i>

```
# Set working directory
cd /scratch/odl7/sturnira_mito

# Set reference paths
REF="references/S_ludovici_QCAZ_18312.fasta"
REFGB="references/S_ludovici_QCAZ_18312.gb"

# Output path
OUTDIR="results"
OUTBED="$OUTDIR/cds_coords.bed"

```

