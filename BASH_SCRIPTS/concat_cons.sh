
#!/usr/bin/env bash

set -euo pipefail

############################################
# PURPOSE
############################################
# Build optional, UNALIGNED complete-mitochondrial-genome multi-FASTA
# matrices for the <=40%, <=30%, and <=20% whole-consensus
# missing-data thresholds.
#
# This is a complementary full-mitogenome analysis. It is not independent
# evidence from the 13-PCG analysis because the complete mitogenome contains
# those same protein-coding genes.
#
# IMPORTANT:
# - This script combines taxa into multi-FASTA files.
# - It does not align sequences.
# - The revised mafft_cons.slurm must align the resulting matrices.
# - Ingroup consensuses must preserve reference coordinates and reference
#   length. Do not use the older indel-shifted consensus files.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
RESULTS_DIR="results"
REF="references/S_ludovici_QCAZ_18312.fasta"

# Expected output location of the revised make_masked_consensus.slurm.
# Keep these names synchronized when that upstream script is revised.
CONSENSUS_SUBDIR="consensus_masked"
CONSENSUS_SUFFIX="_mt_masked.fasta"

# These lists will be produced by the revised filter_from_con.sh.
KEEP40="$RESULTS_DIR/keep_samples_cons_le40.txt"
KEEP30="$RESULTS_DIR/keep_samples_cons_le30.txt"
KEEP20="$RESULTS_DIR/keep_samples_cons_le20.txt"

OUTROOT="$RESULTS_DIR/consensus_matrices"
REPORT="$OUTROOT/concat_cons_validation.tsv"

# Outgroups must:
# 1. contain exactly one complete mitochondrial sequence;
# 2. use the sample name below as the FASTA identifier;
# 3. be in the same strand orientation as the reference;
# 4. be linearized at the same homologous mitochondrial origin.
#
# The revised outgroup_from_gb.sh will be responsible for this contract.
OUTGROUPS=(
    "Artibeus_PP853570.1"
    "Glossophaga_NC_065682.1"
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

REF=$(realpath "$REF")
RESULTS_DIR=$(realpath "$RESULTS_DIR")
KEEP40=$(realpath "$KEEP40")
KEEP30=$(realpath "$KEEP30")
KEEP20=$(realpath "$KEEP20")
OUTROOT=$(realpath -m "$OUTROOT")
REPORT=$(realpath -m "$REPORT")

[[ -s "$REF" ]] || {
    echo "ERROR: reference FASTA not found or empty: $REF" >&2
    exit 1
}

for file in "$KEEP40" "$KEEP30" "$KEEP20"; do
    [[ -s "$file" ]] || {
        echo "ERROR: keep list not found or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$(dirname "$OUTROOT")"

OUTGROUPS_TEXT=$(IFS=,; echo "${OUTGROUPS[*]}")

############################################
# VALIDATE AND BUILD MATRICES
############################################
python3 - \
    "$REF" \
    "$RESULTS_DIR" \
    "$CONSENSUS_SUBDIR" \
    "$CONSENSUS_SUFFIX" \
    "$KEEP40" \
    "$KEEP30" \
    "$KEEP20" \
    "$OUTROOT" \
    "$REPORT" \
    "$OUTGROUPS_TEXT" <<'PY'
import csv
import shutil
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


reference_path = Path(sys.argv[1])
results_dir = Path(sys.argv[2])
consensus_subdir = sys.argv[3]
consensus_suffix = sys.argv[4]
keep_paths = {
    40: Path(sys.argv[5]),
    30: Path(sys.argv[6]),
    20: Path(sys.argv[7]),
}
outroot = Path(sys.argv[8])
report_path = Path(sys.argv[9])
outgroups = [x for x in sys.argv[10].split(",") if x]

# Ungapped IUPAC DNA symbols. Gaps belong only in the later aligned matrices.
allowed_bases = set("ACGTRYSWKMBDHVN")
unambiguous = set("ACGT")


def read_single_fasta(path):
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one FASTA record; found {}."
            .format(path, len(records))
        )
    return records[0]


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


def load_consensus(sample, is_outgroup):
    path = (
        results_dir
        / sample
        / consensus_subdir
        / (sample + consensus_suffix)
    )

    if not path.is_file():
        raise FileNotFoundError(
            "Missing consensus for {}: {}".format(sample, path)
        )

    record = read_single_fasta(path)

    if record.id != sample:
        raise ValueError(
            "FASTA identifier '{}' in {} does not exactly match expected "
            "sample name '{}'."
            .format(record.id, path, sample)
        )

    sequence = str(record.seq).upper().replace("U", "T")

    if not sequence:
        raise ValueError("Empty sequence in {}.".format(path))

    invalid = sorted(set(sequence) - allowed_bases)
    if invalid:
        raise ValueError(
            "Invalid or aligned symbols in {}: {}. Consensus inputs must be "
            "ungapped IUPAC DNA."
            .format(path, ",".join(invalid))
        )

    if not is_outgroup and len(sequence) != reference_length:
        raise ValueError(
            "Ingroup {} has length {} bp; expected the fixed reference "
            "length of {} bp. This suggests an old indel-shifted consensus "
            "or an upstream coordinate error."
            .format(sample, len(sequence), reference_length)
        )

    n_count = sequence.count("N")
    ambiguity_count = sum(
        1 for base in sequence
        if base not in unambiguous and base != "N"
    )
    called_count = len(sequence) - n_count
    gc_count = sequence.count("G") + sequence.count("C")
    gc_percent_called = (
        100.0 * gc_count / called_count if called_count else 0.0
    )

    stats = {
        "taxon": sample,
        "type": "outgroup" if is_outgroup else "ingroup",
        "length_nt": len(sequence),
        "N_count": n_count,
        "percent_N": 100.0 * n_count / len(sequence),
        "other_ambiguity_count": ambiguity_count,
        "GC_percent_among_non_N": gc_percent_called,
        "source_file": str(path),
    }

    return sequence, stats


reference_record = read_single_fasta(reference_path)
reference_sequence = str(reference_record.seq).upper().replace("U", "T")
reference_length = len(reference_sequence)

if not reference_sequence:
    raise ValueError("Reference FASTA sequence is empty.")

reference_invalid = sorted(set(reference_sequence) - allowed_bases)
if reference_invalid:
    raise ValueError(
        "Reference contains invalid symbols: {}."
        .format(",".join(reference_invalid))
    )

keep = {
    threshold: read_keep_list(path)
    for threshold, path in keep_paths.items()
}

# The sensitivity datasets must be nested.
if not set(keep[20]).issubset(set(keep[30])):
    raise ValueError(
        "The <=20% whole-consensus keep list is not a subset of the <=30% list."
    )

if not set(keep[30]).issubset(set(keep[40])):
    raise ValueError(
        "The <=30% whole-consensus keep list is not a subset of the <=40% list."
    )

if len(outgroups) != len(set(outgroups)):
    raise ValueError("The OUTGROUPS array contains duplicate names.")

all_ingroup = []
seen_ingroup = set()

for threshold in (40, 30, 20):
    for sample in keep[threshold]:
        if sample not in seen_ingroup:
            seen_ingroup.add(sample)
            all_ingroup.append(sample)

ingroup_outgroup_overlap = sorted(set(all_ingroup) & set(outgroups))
if ingroup_outgroup_overlap:
    raise ValueError(
        "Outgroup names also occur in ingroup keep lists: {}."
        .format(", ".join(ingroup_outgroup_overlap))
    )

# Load and validate each taxon once.
sequence_by_taxon = {}
stats_by_taxon = {}

for sample in all_ingroup:
    sequence, stats = load_consensus(sample, is_outgroup=False)
    sequence_by_taxon[sample] = sequence
    stats_by_taxon[sample] = stats

for sample in outgroups:
    sequence, stats = load_consensus(sample, is_outgroup=True)
    sequence_by_taxon[sample] = sequence
    stats_by_taxon[sample] = stats

# Confirm that each keep list still agrees with current consensus missingness.
# A small floating-point tolerance avoids boundary artifacts.
tolerance = 1e-9

for threshold in (40, 30, 20):
    for sample in keep[threshold]:
        percent_n = stats_by_taxon[sample]["percent_N"]
        if percent_n > threshold + tolerance:
            raise ValueError(
                "{} is listed in the <={}% dataset but currently has "
                "{:.6f}% N."
                .format(sample, threshold, percent_n)
            )

# Write everything to a temporary directory first.
temp_parent = outroot.parent
temp_dir = Path(tempfile.mkdtemp(
    prefix=outroot.name + ".tmp.",
    dir=str(temp_parent)
))

try:
    report_rows = [
        stats_by_taxon[taxon]
        for taxon in all_ingroup + outgroups
    ]

    with (temp_dir / report_path.name).open("w", newline="") as handle:
        fieldnames = [
            "taxon",
            "type",
            "length_nt",
            "N_count",
            "percent_N",
            "other_ambiguity_count",
            "GC_percent_among_non_N",
            "source_file",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(report_rows)

    for threshold in (40, 30, 20):
        label = "M{}".format(threshold)
        matrix_dir = temp_dir / label
        matrix_dir.mkdir(parents=True, exist_ok=True)

        taxa = keep[threshold] + outgroups
        output_matrix = (
            matrix_dir
            / "Sturnira_complete_mt_{}.unaligned.fasta".format(label)
        )

        records = [
            SeqRecord(
                Seq(sequence_by_taxon[taxon]),
                id=taxon,
                description=""
            )
            for taxon in taxa
        ]

        written_count = SeqIO.write(records, str(output_matrix), "fasta")
        if written_count != len(taxa):
            raise ValueError(
                "Wrote {} records to {}; expected {}."
                .format(written_count, output_matrix, len(taxa))
            )

        reread = list(SeqIO.parse(str(output_matrix), "fasta"))
        if [record.id for record in reread] != taxa:
            raise ValueError(
                "Taxon order changed while writing {}."
                .format(output_matrix)
            )

        (matrix_dir / "taxa.txt").write_text(
            "".join(taxon + "\n" for taxon in taxa)
        )

        ingroup_lengths = {
            len(sequence_by_taxon[taxon])
            for taxon in keep[threshold]
        }
        outgroup_lengths = [
            len(sequence_by_taxon[taxon])
            for taxon in outgroups
        ]

        summary_path = matrix_dir / "matrix_summary.tsv"
        with summary_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([
                "threshold_percent_N",
                "reference_length_nt",
                "ingroup_taxa",
                "outgroup_taxa",
                "total_taxa",
                "ingroup_unique_lengths",
                "outgroup_lengths",
                "alignment_status",
            ])
            writer.writerow([
                threshold,
                reference_length,
                len(keep[threshold]),
                len(outgroups),
                len(taxa),
                ",".join(map(str, sorted(ingroup_lengths))),
                ",".join(map(str, outgroup_lengths)),
                "unaligned",
            ])

    # Replace old output only after all matrices pass.
    if outroot.exists():
        shutil.rmtree(str(outroot))

    shutil.move(str(temp_dir), str(outroot))

except Exception:
    shutil.rmtree(str(temp_dir), ignore_errors=True)
    raise

print("Reference ID: {}".format(reference_record.id))
print("Reference length: {} bp".format(reference_length))
print("Validated ingroup taxa: {}".format(len(all_ingroup)))
print("Validated outgroups: {}".format(len(outgroups)))
print("M40 ingroup taxa: {}".format(len(keep[40])))
print("M30 ingroup taxa: {}".format(len(keep[30])))
print("M20 ingroup taxa: {}".format(len(keep[20])))
print("Output root: {}".format(outroot))
print("Validation report: {}".format(report_path))
PY

############################################
# FINAL SHELL-LEVEL CHECKS
############################################
for threshold in 40 30 20; do
    LABEL="M${threshold}"
    MATRIX="$OUTROOT/$LABEL/Sturnira_complete_mt_${LABEL}.unaligned.fasta"
    TAXA_FILE="$OUTROOT/$LABEL/taxa.txt"
    SUMMARY="$OUTROOT/$LABEL/matrix_summary.tsv"

    for file in "$MATRIX" "$TAXA_FILE" "$SUMMARY"; do
        [[ -s "$file" ]] || {
            echo "ERROR: expected output missing or empty: $file" >&2
            exit 1
        }
    done

    N_TAXA=$(grep -cv '^[[:space:]]*$' "$TAXA_FILE")
    N_SEQUENCES=$(grep -c '^>' "$MATRIX")

    if [[ "$N_SEQUENCES" -ne "$N_TAXA" ]]; then
        echo "ERROR: $MATRIX has $N_SEQUENCES sequences; expected $N_TAXA." >&2
        exit 1
    fi

    echo "$LABEL: $N_SEQUENCES complete mitochondrial sequences, unaligned."
done

[[ -s "$REPORT" ]] || {
    echo "ERROR: validation report missing or empty: $REPORT" >&2
    exit 1
}

echo
echo "Whole-mitogenome matrix construction completed successfully."
echo "Next step: inspect outgroup orientation/origin and run the revised mafft_cons.slurm."
