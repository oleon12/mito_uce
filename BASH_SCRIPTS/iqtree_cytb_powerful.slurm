#!/usr/bin/env bash
#SBATCH --job-name=cytb_iqtree_Sturnira
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/cytb_iqtree_cytb_iqtree_Sturnira%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/cytb_iqtree_cytb_iqtree_Sturnira%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

# Publication-oriented IQ-TREE analysis for the Sturnira cyt-b alignment.
# Tailored to Sturnira_cytb.mafft.fasta inspected on 2026-07-31.
#
# The uploaded alignment has 1152 columns. Cyt-b occupies columns 1-1140;
# columns 1141-1152 are downstream flanking sequence in the outgroups.
# This script trims those columns, performs QC, creates codon partitions,
# tests symmetry, and runs a thorough partitioned maximum-likelihood analysis.

set -euo pipefail

WORKDIR="/scratch/odl7/sturnira_mito"
RAW_ALIGN="$WORKDIR/results/alignments/Sturnira_cytb.mafft.fasta"
OUTDIR="$WORKDIR/results/"
CLEAN_ALIGN="$OUTDIR/alignments/Sturnira_cytb.trimmed_1140.fasta"
PARTITIONS="$OUTDIR/iqtree_FINAL/cytb_codon_partitions.txt"
QC_REPORT="$OUTDIR/iqtree_FINAL/Sturnira_cytb_QC.tsv"
OUTGROUP="Artibeus_PP853570.1,Glossophaga_NC_065682.1"
PREFIX="$OUTDIR/iqtree_FINAL/Sturnira_cytb_partitioned"
SEED=20260731
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate iqtree

# Prefer IQ-TREE 3, but remain compatible with IQ-TREE 2 installations.
IQTREE=""
for candidate in iqtree3 iqtree2 iqtree; do
    if command -v "$candidate" >/dev/null 2>&1; then
        IQTREE="$(command -v "$candidate")"
        break
    fi
done

if [[ -z "$IQTREE" ]]; then
    echo "ERROR: IQ-TREE executable not found in conda environment 'iqtree'." >&2
    exit 1
fi

if [[ ! -s "$RAW_ALIGN" ]]; then
    echo "ERROR: alignment not found or empty: $RAW_ALIGN" >&2
    exit 1
fi

echo "IQ-TREE executable: $IQTREE"
"$IQTREE" --version
echo "Input alignment: $RAW_ALIGN"

# ------------------------------------------------------------
# 1. Validate and trim the exact alignment to cyt-b (1140 bp)
# ------------------------------------------------------------
# Uses only the Python standard library; Biopython is not required.
# The embedded code is compatible with both Python 2.7 and Python 3.
PYTHON_BIN=""
for candidate in python3 python; do
    if command -v "$candidate" >/dev/null 2>&1; then
        PYTHON_BIN="$(command -v "$candidate")"
        break
    fi
done

if [[ -z "$PYTHON_BIN" ]]; then
    echo "ERROR: neither python3 nor python is available." >&2
    exit 1
fi

echo "Python executable: $PYTHON_BIN"
"$PYTHON_BIN" --version 2>&1

"$PYTHON_BIN" - "$RAW_ALIGN" "$CLEAN_ALIGN" "$QC_REPORT" <<'PY'
from __future__ import print_function
import sys
from collections import Counter

input_path = sys.argv[1]
output_path = sys.argv[2]
qc_path = sys.argv[3]

def read_fasta(path):
    records = []
    name = None
    chunks = []
    with open(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records

records = read_fasta(input_path)
if not records:
    raise SystemExit("ERROR: no FASTA records detected.")

names = [name for name, seq in records]
if len(names) != len(set(names)):
    duplicates = [name for name, n in Counter(names).items() if n > 1]
    raise SystemExit("ERROR: duplicated sequence names: {0}".format(duplicates))

lengths = sorted(set(len(seq) for name, seq in records))
if len(lengths) != 1:
    raise SystemExit(
        "ERROR: sequences have unequal aligned lengths: {0}".format(lengths)
    )

alignment_length = lengths[0]
if alignment_length != 1152:
    raise SystemExit(
        "ERROR: expected the inspected 1152-column alignment, found {0}. "
        "Recheck coordinates before trimming.".format(alignment_length)
    )

allowed = set("ACGTRYSWKMBDHVN?-")
bad = sorted(set(
    char
    for name, seq in records
    for char in seq
    if char not in allowed
))
if bad:
    raise SystemExit("ERROR: unsupported sequence symbols: {0}".format(bad))

trimmed = [(name, seq[:1140]) for name, seq in records]

with open(output_path, "w") as handle:
    for name, seq in trimmed:
        handle.write(">{0}\n".format(name))
        for start in range(0, len(seq), 60):
            handle.write(seq[start:start + 60] + "\n")

# Vertebrate mitochondrial stop codons.
stop_codons = set(["TAA", "TAG", "AGA", "AGG"])
fieldnames = [
    "taxon",
    "alignment_length",
    "non_gap_sites",
    "coverage_percent",
    "N_gap_question_sites",
    "other_IUPAC_ambiguities",
    "partial_gap_codons",
    "internal_stop_codons",
    "terminal_stop_codon",
]
rows = []
internal_stop_taxa = []

for name, seq in trimmed:
    codons = [seq[i:i + 3] for i in range(0, 1140, 3)]
    internal_stops = [
        str(index + 1)
        for index, codon in enumerate(codons[:-1])
        if codon in stop_codons
    ]
    terminal_stop = codons[-1] if codons[-1] in stop_codons else ""
    if internal_stops:
        internal_stop_taxa.append((name, internal_stops))

    non_gap = sum(1 for base in seq if base != "-")
    row = {
        "taxon": name,
        "alignment_length": len(seq),
        "non_gap_sites": non_gap,
        "coverage_percent": "{0:.2f}".format(
            100.0 * non_gap / float(len(seq))
        ),
        "N_gap_question_sites": sum(1 for base in seq if base in "N?-"),
        "other_IUPAC_ambiguities": sum(
            1 for base in seq if base not in "ACGTN?-"
        ),
        "partial_gap_codons": sum(
            1 for codon in codons if "-" in codon and codon != "---"
        ),
        "internal_stop_codons": ",".join(internal_stops),
        "terminal_stop_codon": terminal_stop,
    }
    rows.append(row)

with open(qc_path, "w") as handle:
    handle.write("\t".join(fieldnames) + "\n")
    for row in rows:
        handle.write(
            "\t".join(str(row[field]) for field in fieldnames) + "\n"
        )

print("Validated {0} taxa.".format(len(trimmed)))
print("Original alignment length: {0}".format(alignment_length))
print("Cyt-b alignment retained: columns 1-1140")
print("Cleaned alignment written to: {0}".format(output_path))
print("QC report written to: {0}".format(qc_path))

if internal_stop_taxa:
    print("WARNING: internal vertebrate-mitochondrial stop codons detected:")
    for name, positions in internal_stop_taxa:
        print("  {0}: codon(s) {1}".format(name, ",".join(positions)))
    print(
        "Verify these sequences for sequencing error, NUMTs, "
        "or incorrect base calls."
    )
PY

# ------------------------------------------------------------
# 2. Create three initial codon-position partitions
# ------------------------------------------------------------
cat > "$PARTITIONS" <<'EOF'
DNA, cytb_pos1 = 1-1140\3
DNA, cytb_pos2 = 2-1140\3
DNA, cytb_pos3 = 3-1140\3
EOF

# ------------------------------------------------------------
# 3. Test stationarity and homogeneity by partition
# ------------------------------------------------------------
"$IQTREE" \
    -s "$CLEAN_ALIGN" \
    -st DNA \
    -p "$PARTITIONS" \
    --symtest-only \
    -T "$SLURM_CPUS_PER_TASK" \
    --prefix "$OUTDIR/Sturnira_cytb_symmetry"

# ------------------------------------------------------------
# 4. Main partitioned ML analysis
# ------------------------------------------------------------
# MFP+MERGE selects the best model for each initial codon partition
# and tests whether partitions should be merged.
#
# -B 5000       : 5,000 ultrafast bootstrap replicates
# --alrt 1000   : 1,000 SH-aLRT replicates
# -bnni         : NNI optimization of UFBoot trees
# -allnni       : more thorough main-tree NNI search
# -nstop 500    : longer search before stopping
# -pers 0.2     : appropriate safeguard with many partial sequences
# -wbt          : retain bootstrap trees
#
# IQ-TREE checkpointing is automatic. Re-submit the same script after an
# interruption to resume; do not add -redo unless you intentionally want
# to overwrite and restart the complete analysis.

"$IQTREE" \
    -s "$CLEAN_ALIGN" \
    -st DNA \
    -p "$PARTITIONS" \
    -m MFP+MERGE \
    -B 5000 \
    --alrt 1000 \
    -bnni \
    -allnni \
    -nstop 500 \
    -pers 0.2 \
    -wbt \
    -T "$SLURM_CPUS_PER_TASK" \
    -seed "$SEED" \
    -o "$OUTGROUP" \
    --prefix "$PREFIX"

echo
echo "IQ-TREE cyt-b analysis completed."
echo "Main tree:       ${PREFIX}.treefile"
echo "Consensus tree:  ${PREFIX}.contree"
echo "Full report:     ${PREFIX}.iqtree"
echo "Best scheme:     ${PREFIX}.best_scheme.nex (name may vary by IQ-TREE version)"
echo "Symmetry report: $OUTDIR/Sturnira_cytb_symmetry.symtest.csv"
echo "QC report:       $QC_REPORT"
