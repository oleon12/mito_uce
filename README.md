# Mito - UCE

This is a brief guide or tutorial on obtaining mitochondrial DNA, whether in the form of complete genomes or UCE genes. For this guide to work, your UCE data must be enriched during sequencing. This guide is a modification of the original [pipeline](https://github.com/Agustol/mtDNA-mitofinder-pipeline), so all credits go to [Agusto Luzuriaga-Neira](https://github.com/Agustol)

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
<br>

## 0. Environments

Before you start the guide, you must build specific environments. In fact, you must create two different environments, [one](https://github.com/oleon12/mito_uce/blob/main/CONFS/environment.yml) for samtools, bwa, beftools, bedtools, and [another](https://github.com/oleon12/mito_uce/blob/main/CONFS/mitofinder_env.yml) for Python 2.7, blast, spades, and mitofinder-compatible packages. You need to use conda to build the two environment (.yml) files in the <b>CONFS</b> folder.

```
# Make sure that you are in the main before running these lines.

conda env create -f CONFS/environment.yml
conda env create -f CONFS/mitofinder_env.yml --solver=classic

```
<br>

## 1. Reference genome

To extract the genome or genes from your data, you need a reference genome. In this case, I was working with species of the genus Sturnira, so I went to GenBank and downloaded the .fasta and .gb files for a complete genome of <b><i>Sturnira ludovici</i></b>. Thus, you will need to find a reference genome and save it into the <b>references</b> folder.

### 1.1. Prepare reference 

The first step is to prepare the reference files from the reference genome. For this step, you will run the script [<b>prep_ref.slurm</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/prep_ref.slurm). This script first builds the BWA index, which creates several files: .bwt, .sa, .amb, .ann, and .pac. Then it will create a FASTA index with samtools, producing a .fai file (tab-separated: chromosome name, length, offset, etc.). For this script, you must make sure that you call the proper environment with conda, set the path to the reference fasta file, and set the proper working directory (main).

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

Now, you must run the script [<b>cds_bed.slurm</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/cds_bed.slurm). This script takes a GenBank (.gb) file and a matching FASTA file, extracts all CDS (coding sequence) coordinates, and converts them into BED format (0‑based, half‑open) for downstream analyses. For this script, you make sure to set the working directory, the fasta file, and the gb file of the reference genome. The output files will be saved in the <b>results</b> folder. <i>You can change the output file name if you want, but you will need to update it in subsequent steps.</i>

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
<br>

## 2. Mapping reference

Now you will run the script [<b>map_all_ref.slurm</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/map_all_ref.slurm). This script maps paired-end reads (FASTQ) to a reference genome using BWA‑MEM, sorts the resulting BAM file, and filters to keep only mapped reads. This script uses the <b>fastq.gz</b> files and also a text file named [<b>sample_list.txt</b>](https://github.com/oleon12/mito_uce/blob/main/CONFS/sample_list.txt). This text file, in the <b>CONFS</b> folder, contains the match fastq.gz file for each species. You need to create it and save it in the specific folder.

Example of <b>sample_list.txt</b> (be aware that the path to the fastq.gz file must be included):

```
S_angeli_AMNH_214197 raw_data/UCE/CBT03_L0074_R1.fastq.gz raw_data/UCE/CBT03_L0074_R2.fastq.gz
S_bogotensis_AMNH_246573 raw_data/UCE/CBT03_L0096_R1.fastq.gz raw_data/UCE/CBT03_L0096_R2.fastq.gz
```

Once you have finished the text file, you go to the <b>map_all_ref.slurm</b> set up the working directory and the required files and paths. All of the results will be saved in the <b>results</b> folder.

```
# Set up the working directory
cd /scratch/odl7/sturnira_mito

# Set up the sample list file and the reference fasta
SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"

```
### 2.1. Test Mapping

Once you have finished the mapping, you can quickly check your results by running the script [<b>test_map_ref.sh</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/test_map_ref.sh). This will give you the mapped reads, covered position, and average depth.

```
# Run in the terminal

sh BASH_SCRIPTS/test_map_ref.sh

```

### 2.2. Mapping summary

Like the previous one, the script [<b>mapping_summary.slurm</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/mapping_summary.slurm) provides a summary of the mapping and saves everything in a TSV file in the <b>results</b> folder. This is an example of a summary table.

```
Sample	                  Mapped_Reads	  Covered_Positions	  Average_Depth
S_angeli_AMNH_213959          420               12318            3.6620
S_angeli_AMNH_214197          293               4873             4.3567
S_bogotensis_AMNH_207853      3900              16286            15.6029
S_bogotensis_AMNH_207854      2975              15947            10.0561
S_bogotensis_AMNH_207855      3901              15499            11.4085
```
<br>

## 3. Make VCF

Now, you run the script [<b>vcf_file.slurm</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/vcf_file.slurm) which processes all mapped samples to produce VCF files. A VCF file (Variant Call Format) is a standard text format used in bioinformatics to store genetic variants (differences) between a sample's genome and a reference genome. Similar to the previous script, you will need the <b>sample_list.txt</b> file and the <b> reference genome</b>.

```
# Set up the working directory
cd /scratch/odl7/sturnira_mito

# Set up the sample list file and the reference fasta

SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"

```
### 3.1. Test VCF

Once you generate all VCF files, you can get the number of variants per sample using the script [<b>test_vcf.sh</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/test_vcf.sh).

```
# Run in the terminal

sh BASH_SCRIPTS/test_vcf.sh
```

---

Once you have completed these steps, you will find in the <b>results</b> folder a specific folder for each sample/species, containing the associated <b>BAM</b> and <b>VCF</b>. Now, you can proceed with the next steps. 

```
results/
  ├── S_bogotensis_AMNH_207854/
  │   ├── bam/
  │   │   ├── S_bogotensis_AMNH_207854.sorted.bam
  │   │   ├── S_bogotensis_AMNH_207854.sorted.bam.bai
  │   │   ├── S_bogotensis_AMNH_207854.mapped.bam
  │   │   └── S_bogotensis_AMNH_207854.mapped.bam.bai
  │   └── vcf/
  │       ├── S_bogotensis_AMNH_207854.vcf.gz
  │       └── S_bogotensis_AMNH_207854.vcf.gz.csi
  ├── S_ludovici_QCAZ_18312/
  │   └── ...
  └── ...

```

---

<br>

## 4. Make Consensus

Now with the <b>VCF</b> and <b>BAM</b> files, you can create a consensus fasta file using the script [<b>make_consensus.slurm</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/make_consensus.slurm). For this script, you will need the <b>sample_list.txt</b> and <b>reference genome</b> files. All results will be saved in the <b>results</b> folder. The script automatically handles the <b>VCF</b> files for every sample/species. If you change any name from previous scripts, you will need to modify the script; if not, you are ready to go.

```
# Set the environment, sample list, and fasta files
CONDA_ENV="mt_pipeline"
SAMPLE_LIST="CONFS/sample_list.txt"
REF="references/S_ludovici_QCAZ_18312.fasta"
OUTDIR="results"

# In this part, the script finds the VCF files using the sample list
# and save every consensus in individual folders
VCF="$OUTDIR/$SAMPLE/vcf/$SAMPLE.vcf.gz"
CONSDIR="$OUTDIR/$SAMPLE/consensus"
CONS="$CONSDIR/${SAMPLE}_mt.fasta"

```

### 4.1. Test consensus

Once you have finished the consensus, you can check your results using the script [<b>test_consensus.sh</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/test_consensus.sh). This script will calculate the length of the sequences, count the number of N's, and the percentage of them. The information will be saved in a <b>TSV</b> file in the <b>results</b> folder.

```
# Run in the terminal
sh BASH_SCRIPTS/test_consensus.sh

cat results/consensus_stats.tsv

Sample	                  FASTA	                                                                        Length	Ns	Percent_N
S_bogotensis_AMNH_207853	results/S_bogotensis_AMNH_207853/consensus/S_bogotensis_AMNH_207853_mt.fasta	16634	  0	  0.0000
S_bogotensis_AMNH_207854	results/S_bogotensis_AMNH_207854/consensus/S_bogotensis_AMNH_207854_mt.fasta	16640	  0	  0.0000

```

## 5. Masked consensus

Now, you can create an advanced consensus using the script [<b>make_masked_consensus.slurm</b>](https://github.com/oleon12/mito_uce/blob/main/BASH_SCRIPTS/make_masked_consensus.slurm). This is an advanced consensus generation step that masks low‑coverage regions (converts them to Ns) and validates output length. Similar to the previous one, you will need the <b>sample_list.txt</b> and <b>reference genome</b> files. The script automatically handles the <b>VCF</b> and <b>BAM</b> files, and all results will be saved in the <b>results</b> folder. The parameter <b>MIN_DEPTH=3</b>

### 5.1. Test masked consensus
