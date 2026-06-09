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

## 1. Reference genome

To extract the genome or genes from your data, you need a reference genome. In this case, I was working with species of the genus Sturnira, so I went to GenBank and downloaded the .fasta and .gb files for a complete genome of <b><i>Sturnira ludovici</i></b>. Thus, you will need to find a reference genome and save it into the <b>references</b> folder.

### 1.1. Prepare reference 

The first step is to prepare the reference files from the reference genome. For this step, you will run the script <b>prep_ref.slurm</b>. This script first builds the BWA index, which creates several files: .bwt, .sa, .amb, .ann, and .pac. Then it will create a FASTA index using samtools, which produces a .fai file (tab-separated: chromosome name, length, offset, etc.). For this script, you must make sure that you call the proper environment with conda, set the path to the reference fasta file, and set the proper working directory (main).

All files will be save into the <b>references</b> folder.

```
# Set working directory
cd /scratch/odl7/sturnira_mito

# Call environment
CONDA_ENV="mt_pipeline"

# Set references path
REF="references/S_ludovici_QCAZ_18312.fasta"

```

### 1.2. Build CDS beds
