
#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Build one UNALIGNED multi-FASTA matrix per mitochondrial protein-coding gene
# for the <=40%, <=30%, and <=20% missing-data taxon sets.
#
# IMPORTANT:
# - This script does NOT concatenate genes.
# - Each gene must be aligned independently before concatenation.
# - The downstream alignment/concatenation script will consume these outputs.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
RESULTS_DIR="results"

# This file will be produced by the revised extract_cds.slurm.
# It must contain exactly 13 records with headers such as:
#   >SAMPLE_ID|ND1
#   >SAMPLE_ID|ND2
#   ...
#   >SAMPLE_ID|CYTB
GENE_SUBDIR="genes_masked"
GENE_SUFFIX="_13PCG_masked.fasta"

# These lists will be produced by the revised filter_from_cds.sh.
KEEP40="$RESULTS_DIR/keep_samples_cds_le40.txt"
KEEP30="$RESULTS_DIR/keep_samples_cds_le30.txt"
KEEP20="$RESULTS_DIR/keep_samples_cds_le20.txt"

OUTROOT="$RESULTS_DIR/cds_gene_matrices"
REPORT="$OUTROOT/concat_cds_validation.tsv"

# Outgroups are included in every threshold matrix and are not filtered by
# ingroup missing-data thresholds. Their gene files must follow the same format.
OUTGROUPS=(
    "Artibeus_PP853570.1"
    "Glossophaga_NC_065682.1"
)

EXPECTED_GENES=(
    "ND1"
    "ND2"
    "COX1"
    "COX2"
    "ATP8"
    "ATP6"
    "COX3"
    "ND3"
    "ND4L"
    "ND4"
    "ND5"
    "ND6"
    "CYTB"
)

############################################
# SETUP
############################################
cd "$WORKDIR"

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 was not found." >&2
    exit 1
}

python3 - <<'PY_CHECK'
import sys

if sys.version_info < (3, 6):
    raise SystemExit("ERROR: Python 3.6 or newer is required.")

try:
    import Bio
except ImportError as exc:
    raise SystemExit(
        "ERROR: Biopython is required. Install it in the active environment "
        "with: conda install -c conda-forge biopython"
    ) from exc
PY_CHECK

RESULTS_DIR=$(realpath "$RESULTS_DIR")
KEEP40=$(realpath "$KEEP40")
KEEP30=$(realpath "$KEEP30")
KEEP20=$(realpath "$KEEP20")
OUTROOT=$(realpath -m "$OUTROOT")
REPORT=$(realpath -m "$REPORT")

for file in "$KEEP40" "$KEEP30" "$KEEP20"; do
    [[ -s "$file" ]] || {
        echo "ERROR: keep list not found or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$OUTROOT"

OUTGROUPS_TEXT=$(IFS=,; echo "${OUTGROUPS[*]}")
GENES_TEXT=$(IFS=,; echo "${EXPECTED_GENES[*]}")

############################################
# VALIDATE INPUTS AND BUILD GENE MATRICES
############################################
python3 - \
    "$RESULTS_DIR" \
    "$GENE_SUBDIR" \
    "$GENE_SUFFIX" \
    "$KEEP40" \
    "$KEEP30" \
    "$KEEP20" \
    "$OUTROOT" \
    "$REPORT" \
    "$OUTGROUPS_TEXT" \
    "$GENES_TEXT" <<'PY'
import csv
import shutil
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


results_dir = Path(sys.argv[1])
gene_subdir = sys.argv[2]
gene_suffix = sys.argv[3]
keep_paths = {
    40: Path(sys.argv[4]),
    30: Path(sys.argv[5]),
    20: Path(sys.argv[6]),
}
outroot = Path(sys.argv[7])
report_path = Path(sys.argv[8])
outgroups = [x for x in sys.argv[9].split(",") if x]
expected_genes = [x for x in sys.argv[10].split(",") if x]

allowed_bases = set("ACGTRYSWKMBDHVN")


def read_keep_list(path):
    taxa = []
    seen = set()

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            sample = line.split()[0]

            if sample in seen:
                raise ValueError(
                    "Duplicate sample '{}' in {} at line {}."
                    .format(sample, path, line_number)
                )

            seen.add(sample)
            taxa.append(sample)

    if not taxa:
        raise ValueError("No samples were found in {}.".format(path))

    return taxa


def gene_from_record_id(record_id, sample):
    """
    Accept either:
      >SAMPLE|GENE
    or:
      >GENE

    The publication pipeline will use >SAMPLE|GENE.
    """
    if "|" in record_id:
        prefix, gene = record_id.rsplit("|", 1)
        if prefix != sample:
            raise ValueError(
                "Record '{}' does not match expected sample '{}'."
                .format(record_id, sample)
            )
        return gene.upper()

    return record_id.upper()


def load_sample_genes(sample):
    path = results_dir / sample / gene_subdir / (sample + gene_suffix)

    if not path.is_file():
        raise FileNotFoundError(
            "Missing 13-PCG FASTA for {}: {}".format(sample, path)
        )

    records = list(SeqIO.parse(str(path), "fasta"))

    if not records:
        raise ValueError("No FASTA records found in {}.".format(path))

    genes = {}

    for record in records:
        gene = gene_from_record_id(record.id, sample)

        if gene not in expected_genes:
            raise ValueError(
                "Unexpected gene '{}' in {}. Expected only: {}."
                .format(gene, path, ", ".join(expected_genes))
            )

        if gene in genes:
            raise ValueError(
                "Duplicate gene '{}' in {}.".format(gene, path)
            )

        sequence = str(record.seq).upper().replace("U", "T")

        if not sequence:
            raise ValueError(
                "Empty sequence for {}|{} in {}."
                .format(sample, gene, path)
            )

        invalid = sorted(set(sequence) - allowed_bases)
        if invalid:
            raise ValueError(
                "Invalid nucleotide symbols for {}|{} in {}: {}."
                .format(sample, gene, path, ",".join(invalid))
            )

        if len(sequence) % 3 != 0:
            raise ValueError(
                "{}|{} has length {}, which is not divisible by 3."
                .format(sample, gene, len(sequence))
            )

        protein = str(Seq(sequence).translate(table=2))
        internal_stops = protein.count("*")

        if internal_stops:
            raise ValueError(
                "{}|{} contains {} internal stop codon(s)."
                .format(sample, gene, internal_stops)
            )

        genes[gene] = sequence

    missing = [gene for gene in expected_genes if gene not in genes]
    if missing:
        raise ValueError(
            "{} is missing genes: {}."
            .format(sample, ", ".join(missing))
        )

    if len(genes) != len(expected_genes):
        raise ValueError(
            "{} contains {} genes; expected {}."
            .format(sample, len(genes), len(expected_genes))
        )

    return path, genes


keep = {threshold: read_keep_list(path)
        for threshold, path in keep_paths.items()}

# These should be nested sensitivity datasets.
if not set(keep[20]).issubset(set(keep[30])):
    raise ValueError(
        "The <=20% keep list is not a subset of the <=30% keep list."
    )

if not set(keep[30]).issubset(set(keep[40])):
    raise ValueError(
        "The <=30% keep list is not a subset of the <=40% keep list."
    )

all_ingroup = []
seen_ingroup = set()
for threshold in (40, 30, 20):
    for sample in keep[threshold]:
        if sample not in seen_ingroup:
            seen_ingroup.add(sample)
            all_ingroup.append(sample)

overlap = sorted(set(all_ingroup) & set(outgroups))
if overlap:
    raise ValueError(
        "Outgroup names also occur in the ingroup keep lists: {}."
        .format(", ".join(overlap))
    )

if len(outgroups) != len(set(outgroups)):
    raise ValueError("The OUTGROUPS array contains duplicate names.")

# Validate every taxon once before writing any matrix.
all_taxa = all_ingroup + outgroups
sample_data = {}
validation_rows = []

for sample in all_taxa:
    source_path, genes = load_sample_genes(sample)
    sample_data[sample] = genes
    sample_type = "outgroup" if sample in outgroups else "ingroup"

    for gene in expected_genes:
        sequence = genes[gene]
        n_count = sequence.count("N")
        validation_rows.append({
            "taxon": sample,
            "type": sample_type,
            "gene": gene,
            "length_nt": len(sequence),
            "N_count": n_count,
            "percent_N": 100.0 * n_count / len(sequence),
            "internal_stop_count": 0,
            "source_file": str(source_path),
        })

# Write to a temporary directory first. Incomplete matrices are never retained.
temp_parent = outroot.parent
temp_dir = Path(tempfile.mkdtemp(
    prefix=outroot.name + ".tmp.",
    dir=str(temp_parent)
))

try:
    for threshold in (40, 30, 20):
        label = "M{}".format(threshold)
        matrix_dir = temp_dir / label / "unaligned_by_gene"
        matrix_dir.mkdir(parents=True, exist_ok=True)

        taxa = keep[threshold] + outgroups

        (temp_dir / label / "taxa.txt").write_text(
            "".join(taxon + "\n" for taxon in taxa)
        )

        for gene in expected_genes:
            records = [
                SeqRecord(
                    Seq(sample_data[taxon][gene]),
                    id=taxon,
                    description=""
                )
                for taxon in taxa
            ]

            output = matrix_dir / (gene + ".fasta")
            SeqIO.write(records, str(output), "fasta")

            written = list(SeqIO.parse(str(output), "fasta"))
            if len(written) != len(taxa):
                raise ValueError(
                    "{} contains {} records; expected {}."
                    .format(output, len(written), len(taxa))
                )

            if [record.id for record in written] != taxa:
                raise ValueError(
                    "Taxon order changed while writing {}.".format(output)
                )

        summary_path = temp_dir / label / "matrix_summary.tsv"
        with summary_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([
                "threshold_percent_N",
                "ingroup_taxa",
                "outgroup_taxa",
                "total_taxa",
                "genes",
            ])
            writer.writerow([
                threshold,
                len(keep[threshold]),
                len(outgroups),
                len(taxa),
                len(expected_genes),
            ])

    report_temp = temp_dir / report_path.name
    with report_temp.open("w", newline="") as handle:
        fieldnames = [
            "taxon",
            "type",
            "gene",
            "length_nt",
            "N_count",
            "percent_N",
            "internal_stop_count",
            "source_file",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(validation_rows)

    # Replace the old output only after every validation has succeeded.
    if outroot.exists():
        shutil.rmtree(str(outroot))

    shutil.move(str(temp_dir), str(outroot))

except Exception:
    shutil.rmtree(str(temp_dir), ignore_errors=True)
    raise

print("Validated taxa: {}".format(len(all_taxa)))
print("Outgroups: {}".format(", ".join(outgroups)))
print("M40 ingroup taxa: {}".format(len(keep[40])))
print("M30 ingroup taxa: {}".format(len(keep[30])))
print("M20 ingroup taxa: {}".format(len(keep[20])))
print("Genes per matrix: {}".format(len(expected_genes)))
print("Output root: {}".format(outroot))
print("Validation report: {}".format(report_path))
PY

############################################
# FINAL SHELL-LEVEL CHECKS
############################################
for threshold in 40 30 20; do
    MATRIX_DIR="$OUTROOT/M${threshold}/unaligned_by_gene"
    TAXA_FILE="$OUTROOT/M${threshold}/taxa.txt"

    [[ -d "$MATRIX_DIR" ]] || {
        echo "ERROR: matrix directory missing: $MATRIX_DIR" >&2
        exit 1
    }

    [[ -s "$TAXA_FILE" ]] || {
        echo "ERROR: taxa file missing or empty: $TAXA_FILE" >&2
        exit 1
    }

    N_TAXA=$(grep -cv '^[[:space:]]*$' "$TAXA_FILE")
    N_GENES=$(find "$MATRIX_DIR" -maxdepth 1 -type f -name '*.fasta' | wc -l)

    if [[ "$N_GENES" -ne 13 ]]; then
        echo "ERROR: M${threshold} has $N_GENES gene matrices; expected 13." >&2
        exit 1
    fi

    for gene in "${EXPECTED_GENES[@]}"; do
        FASTA="$MATRIX_DIR/${gene}.fasta"

        [[ -s "$FASTA" ]] || {
            echo "ERROR: missing gene matrix: $FASTA" >&2
            exit 1
        }

        N_RECORDS=$(grep -c '^>' "$FASTA")

        if [[ "$N_RECORDS" -ne "$N_TAXA" ]]; then
            echo "ERROR: $FASTA has $N_RECORDS sequences; expected $N_TAXA." >&2
            exit 1
        fi
    done

    echo "M${threshold}: $N_TAXA taxa across 13 unaligned gene matrices."
done

echo
echo "Gene-matrix construction completed successfully."
echo "No genes were concatenated."
echo "Next step: align each gene independently, then concatenate the alignments."
