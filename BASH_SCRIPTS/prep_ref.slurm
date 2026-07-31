#!/usr/bin/env bash
#SBATCH --job-name=prep_mt_ref
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/odl7/sturnira_mito/LOG/prep_mt_ref_%j.out
#SBATCH --error=/scratch/odl7/sturnira_mito/ERR_OUT/prep_mt_ref_%j.err
#SBATCH --mail-user=odl7@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

############################################
# PURPOSE
############################################
# Validate the mitochondrial reference and its GenBank annotation, then build
# and verify the classic BWA and SAMtools FASTA indexes used by the pipeline.
#
# This script is the reference gatekeeper. It verifies:
# - exactly one FASTA and one GenBank record;
# - exact FASTA/GenBank sequence identity;
# - an unambiguous circular DNA reference;
# - conventional mammalian circular origin at tRNA-Phe position 1;
# - all 13 expected protein-coding genes;
# - codon_start=1, translation table 2, and annotated translations;
# - complete and current BWA index components;
# - a current, functional .fai index;
# - a BWA-MEM smoke-test alignment against the prepared reference.

############################################
# CONFIGURATION
############################################
WORKDIR="/scratch/odl7/sturnira_mito"
CONDA_ENV="mt_pipeline"

REF="references/S_ludovici_QCAZ_18312.fasta"
REFGB="references/S_ludovici_QCAZ_18312.gb"

QC_DIR="results/reference_preparation"
REFERENCE_QC="$QC_DIR/reference_qc.tsv"
REFERENCE_CDS_QC="$QC_DIR/reference_cds_qc.tsv"
INDEX_MANIFEST="$QC_DIR/reference_index_manifest.tsv"
INDEXED_REFERENCE_SHA256="$QC_DIR/indexed_reference.sha256"

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

BWA_INDEX_SUFFIXES=(amb ann bwt pac sa)

############################################
# ENVIRONMENT
############################################
cd "$WORKDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

for command_name in python3 bwa samtools sha256sum; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: $command_name was not found in $CONDA_ENV." >&2
        exit 1
    }
done

python3 - <<'PY_CHECK'
import sys

if sys.version_info < (3, 6):
    raise SystemExit("ERROR: Python 3.6 or newer is required.")

try:
    import Bio
except ImportError as exc:
    raise SystemExit(
        "ERROR: Biopython is required. Install it with: "
        "conda install -c conda-forge biopython"
    ) from exc
PY_CHECK

BWA_VERSION="$(bwa 2>&1 | awk '/Version:/ {print $2; exit}' || true)"
SAMTOOLS_VERSION="$(samtools --version | awk 'NR==1 {print $2; exit}')"

if [[ -z "$BWA_VERSION" ]]; then
    BWA_VERSION="unknown"
fi

if [[ -z "$SAMTOOLS_VERSION" ]]; then
    SAMTOOLS_VERSION="unknown"
fi

echo "BWA version:      $BWA_VERSION"
echo "SAMtools version: $SAMTOOLS_VERSION"

############################################
# RESOLVE AND CHECK PATHS
############################################
REF=$(realpath -m "$REF")
REFGB=$(realpath -m "$REFGB")
QC_DIR=$(realpath -m "$QC_DIR")
REFERENCE_QC=$(realpath -m "$REFERENCE_QC")
REFERENCE_CDS_QC=$(realpath -m "$REFERENCE_CDS_QC")
INDEX_MANIFEST=$(realpath -m "$INDEX_MANIFEST")
INDEXED_REFERENCE_SHA256=$(realpath -m "$INDEXED_REFERENCE_SHA256")

for file in "$REF" "$REFGB"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required reference input missing or empty: $file" >&2
        exit 1
    }
done

mkdir -p "$QC_DIR"

GENES_TEXT=$(IFS=,; echo "${EXPECTED_GENES[*]}")

############################################
# VALIDATE FASTA, GENBANK, AND CDS ANNOTATION
############################################
TMPROOT=$(mktemp -d "$QC_DIR/.prep_ref.XXXXXX")
trap 'rm -rf "$TMPROOT"' EXIT

META_ENV="$TMPROOT/reference_metadata.env"
TMP_REFERENCE_QC="$TMPROOT/reference_qc.tsv"
TMP_REFERENCE_CDS_QC="$TMPROOT/reference_cds_qc.tsv"

python3 - \
    "$REF" \
    "$REFGB" \
    "$TMP_REFERENCE_QC" \
    "$TMP_REFERENCE_CDS_QC" \
    "$META_ENV" \
    "$GENES_TEXT" <<'PY'
import csv
import hashlib
import re
import shlex
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition


fasta_path = Path(sys.argv[1])
genbank_path = Path(sys.argv[2])
reference_qc_path = Path(sys.argv[3])
reference_cds_qc_path = Path(sys.argv[4])
metadata_env_path = Path(sys.argv[5])
expected_genes = [x for x in sys.argv[6].split(",") if x]

mitochondrial_table = CodonTable.unambiguous_dna_by_id[2]
mitochondrial_start_codons = set(mitochondrial_table.start_codons)

GENE_ALIASES = {
    "COI": "COX1",
    "CO1": "COX1",
    "COXI": "COX1",
    "COII": "COX2",
    "CO2": "COX2",
    "COXII": "COX2",
    "COIII": "COX3",
    "CO3": "COX3",
    "COXIII": "COX3",
    "CYB": "CYTB",
    "COB": "CYTB",
    "ATPASE6": "ATP6",
    "ATPASE8": "ATP8",
    "NAD1": "ND1",
    "NAD2": "ND2",
    "NAD3": "ND3",
    "NAD4": "ND4",
    "NAD4L": "ND4L",
    "NAD5": "ND5",
    "NAD6": "ND6",
}


def canonical_gene_name(value):
    normalized = re.sub(r"[^A-Za-z0-9]", "", value.upper())

    if normalized in expected_genes:
        return normalized

    return GENE_ALIASES.get(normalized, normalized)


def read_single(path, fmt):
    records = list(SeqIO.parse(str(path), fmt))

    if len(records) != 1:
        raise ValueError(
            "{} must contain exactly one {} record; found {}."
            .format(path, fmt, len(records))
        )

    return records[0]


def normalize(sequence):
    return str(sequence).upper().replace("U", "T")


def sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def exact_location(feature, gene):
    if feature.location.strand not in (1, -1):
        raise ValueError("{} has an undefined CDS strand.".format(gene))

    parts = list(feature.location.parts)

    if not parts:
        raise ValueError("{} has no CDS location parts.".format(gene))

    for part in parts:
        if not isinstance(part.start, ExactPosition):
            raise ValueError(
                "{} has a fuzzy CDS start: {}."
                .format(gene, part.start)
            )

        if not isinstance(part.end, ExactPosition):
            raise ValueError(
                "{} has a fuzzy CDS end: {}."
                .format(gene, part.end)
            )


def validate_translation(coding_sequence, annotated_translation, gene):
    calculated = str(Seq(coding_sequence).translate(table=2))

    if len(calculated) != len(annotated_translation):
        raise ValueError(
            "{} translated length {} differs from annotation length {}."
            .format(gene, len(calculated), len(annotated_translation))
        )

    for index, (observed, expected) in enumerate(
        zip(calculated, annotated_translation)
    ):
        if observed == expected:
            continue

        if (
            index == 0
            and expected == "M"
            and coding_sequence[:3] in mitochondrial_start_codons
        ):
            continue

        raise ValueError(
            "{} translation mismatch at amino-acid position {}: "
            "calculated {}, annotated {}."
            .format(gene, index + 1, observed, expected)
        )

    if "*" in calculated:
        positions = [
            str(index + 1)
            for index, amino_acid in enumerate(calculated)
            if amino_acid == "*"
        ]
        raise ValueError(
            "{} contains internal stop codon(s) at positions {}."
            .format(gene, ",".join(positions))
        )


fasta_record = read_single(fasta_path, "fasta")
genbank_record = read_single(genbank_path, "genbank")

fasta_sequence = normalize(fasta_record.seq)
genbank_sequence = normalize(genbank_record.seq)

if not re.fullmatch(r"[A-Za-z0-9._-]+", fasta_record.id):
    raise ValueError(
        "Unsafe FASTA identifier '{}'. Use only letters, numbers, periods, "
        "underscores, and hyphens."
        .format(fasta_record.id)
    )

if not fasta_sequence:
    raise ValueError("Reference FASTA sequence is empty.")

invalid = sorted(set(fasta_sequence) - set("ACGT"))

if invalid:
    raise ValueError(
        "Reference must be unambiguous A/C/G/T DNA; observed: {}."
        .format(",".join(invalid))
    )

if fasta_sequence != genbank_sequence:
    raise ValueError(
        "Reference FASTA and GenBank sequences are not exactly identical."
    )

if genbank_record.annotations.get("topology", "").lower() != "circular":
    raise ValueError(
        "GenBank topology is not annotated as circular: {!r}."
        .format(genbank_record.annotations.get("topology"))
    )

if genbank_record.annotations.get("molecule_type", "").upper() != "DNA":
    raise ValueError(
        "GenBank molecule_type is not DNA: {!r}."
        .format(genbank_record.annotations.get("molecule_type"))
    )

# Validate the conventional mammalian linearization origin.
trna_phe_features = []

for feature in genbank_record.features:
    if feature.type != "tRNA":
        continue

    values = []
    for qualifier in ("gene", "product", "locus_tag"):
        values.extend(feature.qualifiers.get(qualifier, []))

    normalized_values = {
        re.sub(r"[^A-Za-z0-9]", "", value.upper())
        for value in values
    }

    if "TRNAPHE" in normalized_values:
        trna_phe_features.append(feature)

if len(trna_phe_features) != 1:
    raise ValueError(
        "Expected exactly one tRNA-Phe feature; found {}."
        .format(len(trna_phe_features))
    )

trna_phe = trna_phe_features[0]

if int(trna_phe.location.start) != 0 or trna_phe.location.strand != 1:
    raise ValueError(
        "Reference is not linearized at plus-strand tRNA-Phe position 1: {}."
        .format(trna_phe.location)
    )

features = {}

for feature in genbank_record.features:
    if feature.type != "CDS":
        continue

    candidate_values = []

    for qualifier in ("gene", "locus_tag", "product"):
        candidate_values.extend(feature.qualifiers.get(qualifier, []))

    gene = None

    for value in candidate_values:
        candidate = canonical_gene_name(value)

        if candidate in expected_genes:
            gene = candidate
            break

    if gene is None:
        continue

    exact_location(feature, gene)

    if gene in features:
        raise ValueError(
            "Duplicate reference CDS annotation for {}."
            .format(gene)
        )

    features[gene] = feature

missing = [gene for gene in expected_genes if gene not in features]

if missing:
    raise ValueError(
        "Reference is missing expected CDS annotations: {}."
        .format(", ".join(missing))
    )

cds_rows = []
total_coding_nt = 0

for gene in expected_genes:
    feature = features[gene]
    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    translation_table = int(
        feature.qualifiers.get("transl_table", ["2"])[0]
    )
    annotated_translation = (
        feature.qualifiers.get("translation", [""])[0]
        .replace(" ", "")
        .replace("\n", "")
        .upper()
    )

    if codon_start != 1:
        raise ValueError(
            "{} has codon_start={}; expected 1."
            .format(gene, codon_start)
        )

    if translation_table != 2:
        raise ValueError(
            "{} uses translation table {}; expected 2."
            .format(gene, translation_table)
        )

    if not annotated_translation:
        raise ValueError(
            "{} lacks a /translation qualifier."
            .format(gene)
        )

    full_sequence = normalize(feature.extract(genbank_record.seq))
    coding_nt = len(annotated_translation) * 3
    terminal_stop_nt = len(full_sequence) - coding_nt

    if terminal_stop_nt not in (0, 1, 2, 3):
        raise ValueError(
            "{} full CDS length {} and translation length {} imply an "
            "invalid terminal difference of {} nt."
            .format(
                gene,
                len(full_sequence),
                len(annotated_translation),
                terminal_stop_nt,
            )
        )

    coding_sequence = full_sequence[:coding_nt]

    if len(coding_sequence) % 3 != 0:
        raise ValueError(
            "{} coding length {} is not divisible by 3."
            .format(gene, len(coding_sequence))
        )

    validate_translation(coding_sequence, annotated_translation, gene)
    total_coding_nt += len(coding_sequence)

    cds_rows.append({
        "gene": gene,
        "strand": "+" if feature.location.strand == 1 else "-",
        "location": str(feature.location),
        "full_annotated_nt": len(full_sequence),
        "coding_nt_retained": len(coding_sequence),
        "terminal_stop_nt_excluded": terminal_stop_nt,
        "protein_aa": len(annotated_translation),
        "codon_start": codon_start,
        "translation_table": translation_table,
        "translation_validated": "yes",
    })

reference_sha256 = sha256(fasta_sequence)
gc_percent = 100.0 * (
    fasta_sequence.count("G") + fasta_sequence.count("C")
) / len(fasta_sequence)

with reference_qc_path.open("w", newline="") as handle:
    fieldnames = [
        "fasta_file",
        "genbank_file",
        "fasta_id",
        "genbank_id",
        "organism",
        "topology",
        "molecule_type",
        "length_nt",
        "A_count",
        "C_count",
        "G_count",
        "T_count",
        "N_count",
        "GC_percent",
        "sha256_sequence",
        "fasta_genbank_identical",
        "circular_origin_feature",
        "circular_origin_start_1based",
        "circular_origin_strand",
        "protein_coding_genes",
        "total_coding_nt",
        "reference_validation",
    ]
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerow({
        "fasta_file": str(fasta_path),
        "genbank_file": str(genbank_path),
        "fasta_id": fasta_record.id,
        "genbank_id": genbank_record.id,
        "organism": genbank_record.annotations.get("organism", ""),
        "topology": genbank_record.annotations.get("topology", ""),
        "molecule_type": genbank_record.annotations.get(
            "molecule_type", ""
        ),
        "length_nt": len(fasta_sequence),
        "A_count": fasta_sequence.count("A"),
        "C_count": fasta_sequence.count("C"),
        "G_count": fasta_sequence.count("G"),
        "T_count": fasta_sequence.count("T"),
        "N_count": fasta_sequence.count("N"),
        "GC_percent": gc_percent,
        "sha256_sequence": reference_sha256,
        "fasta_genbank_identical": "yes",
        "circular_origin_feature": "tRNA-Phe",
        "circular_origin_start_1based": 1,
        "circular_origin_strand": "+",
        "protein_coding_genes": len(expected_genes),
        "total_coding_nt": total_coding_nt,
        "reference_validation": "passed",
    })

with reference_cds_qc_path.open("w", newline="") as handle:
    fieldnames = [
        "gene",
        "strand",
        "location",
        "full_annotated_nt",
        "coding_nt_retained",
        "terminal_stop_nt_excluded",
        "protein_aa",
        "codon_start",
        "translation_table",
        "translation_validated",
    ]
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(cds_rows)

metadata = {
    "REF_ID": fasta_record.id,
    "REF_GB_ID": genbank_record.id,
    "REF_LENGTH": str(len(fasta_sequence)),
    "REF_SHA256": reference_sha256,
    "REF_TOTAL_CODING_NT": str(total_coding_nt),
}

with metadata_env_path.open("w") as handle:
    for key, value in metadata.items():
        handle.write("{}={}\n".format(key, shlex.quote(value)))

print("Reference FASTA ID:       {}".format(fasta_record.id))
print("Reference GenBank ID:     {}".format(genbank_record.id))
print("Reference length:         {} bp".format(len(fasta_sequence)))
print("Reference SHA-256:        {}".format(reference_sha256))
print("Protein-coding genes:     {}".format(len(expected_genes)))
print("Retained coding length:   {} bp".format(total_coding_nt))
print("Circular origin:          tRNA-Phe at position 1 (+)")
PY

source "$META_ENV"

CURRENT_FILE_SHA256=$(sha256sum "$REF" | awk '{print $1}')

# Because the FASTA contains one canonical sequence without non-sequence
# characters, the raw-file and parsed-sequence checksums are recorded
# separately. Staleness is tracked using the raw FASTA-file checksum.
PREVIOUS_FILE_SHA256=""

if [[ -s "$INDEXED_REFERENCE_SHA256" ]]; then
    PREVIOUS_FILE_SHA256=$(awk 'NR==1 {print $1}' "$INDEXED_REFERENCE_SHA256")
fi

############################################
# DETERMINE WHETHER INDEXES ARE CURRENT
############################################
REBUILD_BWA="false"
REBUILD_FAI="false"

if [[ "$PREVIOUS_FILE_SHA256" != "$CURRENT_FILE_SHA256" ]]; then
    REBUILD_BWA="true"
    REBUILD_FAI="true"
fi

for suffix in "${BWA_INDEX_SUFFIXES[@]}"; do
    if [[ ! -s "${REF}.${suffix}" ]]; then
        REBUILD_BWA="true"
    fi
done

if [[ ! -s "${REF}.fai" ]]; then
    REBUILD_FAI="true"
else
    FAI_ROWS=$(awk 'NF>0 {n++} END{print n+0}' "${REF}.fai")
    FAI_ID=$(awk 'NR==1 {print $1}' "${REF}.fai")
    FAI_LENGTH=$(awk 'NR==1 {print $2}' "${REF}.fai")

    if [[ "$FAI_ROWS" -ne 1 || "$FAI_ID" != "$REF_ID" || "$FAI_LENGTH" != "$REF_LENGTH" ]]; then
        REBUILD_FAI="true"
    fi
fi

############################################
# BUILD BWA INDEX AT A TEMPORARY PREFIX
############################################
BWA_ACTION="validated_existing"

if [[ "$REBUILD_BWA" == "true" ]]; then
    echo "Building classic BWA index with the IS algorithm..."

    TMP_BWA_PREFIX="$TMPROOT/bwa_reference"
    bwa index -a is -p "$TMP_BWA_PREFIX" "$REF"

    for suffix in "${BWA_INDEX_SUFFIXES[@]}"; do
        [[ -s "${TMP_BWA_PREFIX}.${suffix}" ]] || {
            echo "ERROR: BWA did not create ${TMP_BWA_PREFIX}.${suffix}" >&2
            exit 1
        }

        mv "${TMP_BWA_PREFIX}.${suffix}" "${REF}.${suffix}.new.$$"
    done

    for suffix in "${BWA_INDEX_SUFFIXES[@]}"; do
        mv "${REF}.${suffix}.new.$$" "${REF}.${suffix}"
    done

    BWA_ACTION="rebuilt"
else
    echo "Existing BWA index matches the recorded reference checksum."
fi

############################################
# BUILD SAMTOOLS FAI FROM AN EXACT TEMP COPY
############################################
FAI_ACTION="validated_existing"

if [[ "$REBUILD_FAI" == "true" ]]; then
    echo "Building SAMtools FASTA index..."

    TMP_FASTA="$TMPROOT/reference.fasta"
    cp "$REF" "$TMP_FASTA"
    samtools faidx "$TMP_FASTA"

    [[ -s "${TMP_FASTA}.fai" ]] || {
        echo "ERROR: samtools faidx did not create ${TMP_FASTA}.fai" >&2
        exit 1
    }

    mv "${TMP_FASTA}.fai" "${REF}.fai.new.$$"
    mv "${REF}.fai.new.$$" "${REF}.fai"
    FAI_ACTION="rebuilt"
else
    echo "Existing FASTA index matches the reference ID and length."
fi

############################################
# VERIFY THE FAI BY RETRIEVING THE FULL RECORD
############################################
FETCHED_FASTA="$TMPROOT/reference_from_fai.fasta"
samtools faidx "$REF" "$REF_ID" > "$FETCHED_FASTA"

python3 - "$REF" "$FETCHED_FASTA" <<'PY'
import sys
from Bio import SeqIO

reference = SeqIO.read(sys.argv[1], "fasta")
fetched = SeqIO.read(sys.argv[2], "fasta")

if str(reference.seq).upper() != str(fetched.seq).upper():
    raise SystemExit(
        "ERROR: sequence retrieved through the .fai differs from the reference."
    )

print("SAMtools .fai full-sequence retrieval: passed")
PY

############################################
# VERIFY BWA INDEX WITH AN EXACT REFERENCE READ
############################################
SMOKE_FASTQ="$TMPROOT/reference_index_smoke.fastq"

python3 - "$REF" "$SMOKE_FASTQ" <<'PY'
import sys
from Bio import SeqIO

record = SeqIO.read(sys.argv[1], "fasta")
sequence = str(record.seq[:150]).upper()

if len(sequence) < 70:
    raise SystemExit("ERROR: reference is too short for the BWA-MEM smoke test.")

with open(sys.argv[2], "w") as handle:
    handle.write("@reference_index_smoke\n")
    handle.write(sequence + "\n")
    handle.write("+\n")
    handle.write("I" * len(sequence) + "\n")
PY

SMOKE_MAPPED=$( \
    bwa mem -t 1 "$REF" "$SMOKE_FASTQ" 2> "$TMPROOT/bwa_smoke.log" \
    | samtools view -c -F 4 - \
)

if [[ "$SMOKE_MAPPED" -lt 1 ]]; then
    echo "ERROR: BWA-MEM smoke-test read did not map to the prepared reference." >&2
    cat "$TMPROOT/bwa_smoke.log" >&2
    exit 1
fi

echo "BWA-MEM index smoke test: passed ($SMOKE_MAPPED mapped alignment)"

############################################
# RECORD SUCCESSFUL INDEX STATE AND MANIFEST
############################################
printf '%s  %s\n' "$CURRENT_FILE_SHA256" "$REF" > "$TMPROOT/indexed_reference.sha256"

TMP_INDEX_MANIFEST="$TMPROOT/reference_index_manifest.tsv"

python3 - \
    "$REF" \
    "$REF_ID" \
    "$REF_LENGTH" \
    "$CURRENT_FILE_SHA256" \
    "$BWA_VERSION" \
    "$SAMTOOLS_VERSION" \
    "$BWA_ACTION" \
    "$FAI_ACTION" \
    "$TMP_INDEX_MANIFEST" <<'PY'
import csv
import hashlib
import sys
from pathlib import Path


reference = Path(sys.argv[1])
reference_id = sys.argv[2]
reference_length = int(sys.argv[3])
reference_file_sha256 = sys.argv[4]
bwa_version = sys.argv[5]
samtools_version = sys.argv[6]
bwa_action = sys.argv[7]
fai_action = sys.argv[8]
output = Path(sys.argv[9])

index_files = [
    ("bwa", bwa_version, bwa_action, Path(str(reference) + suffix))
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")
]
index_files.append(
    ("samtools_faidx", samtools_version, fai_action, Path(str(reference) + ".fai"))
)

rows = []

for tool, version, action, path in index_files:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("Missing or empty finalized index: {}.".format(path))

    digest = hashlib.sha256(path.read_bytes()).hexdigest()

    rows.append({
        "reference_file": str(reference),
        "reference_id": reference_id,
        "reference_length_nt": reference_length,
        "reference_file_sha256": reference_file_sha256,
        "tool": tool,
        "tool_version": version,
        "index_action": action,
        "index_file": str(path),
        "index_size_bytes": path.stat().st_size,
        "index_sha256": digest,
        "validation": "passed",
    })

with output.open("w", newline="") as handle:
    fieldnames = [
        "reference_file",
        "reference_id",
        "reference_length_nt",
        "reference_file_sha256",
        "tool",
        "tool_version",
        "index_action",
        "index_file",
        "index_size_bytes",
        "index_sha256",
        "validation",
    ]
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
PY

# Finalize QC and state files only after all validation and smoke tests pass.
mv "$TMP_REFERENCE_QC" "$REFERENCE_QC"
mv "$TMP_REFERENCE_CDS_QC" "$REFERENCE_CDS_QC"
mv "$TMP_INDEX_MANIFEST" "$INDEX_MANIFEST"
mv "$TMPROOT/indexed_reference.sha256" "$INDEXED_REFERENCE_SHA256"

############################################
# FINAL CHECKS
############################################
for file in \
    "$REFERENCE_QC" \
    "$REFERENCE_CDS_QC" \
    "$INDEX_MANIFEST" \
    "$INDEXED_REFERENCE_SHA256" \
    "${REF}.fai"; do

    [[ -s "$file" ]] || {
        echo "ERROR: expected output missing or empty: $file" >&2
        exit 1
    }
done

for suffix in "${BWA_INDEX_SUFFIXES[@]}"; do
    [[ -s "${REF}.${suffix}" ]] || {
        echo "ERROR: BWA index component missing or empty: ${REF}.${suffix}" >&2
        exit 1
    }
done

echo
echo "Reference preparation completed successfully."
echo "Reference ID:        $REF_ID"
echo "Reference length:    $REF_LENGTH bp"
echo "Retained 13-PCG:     $REF_TOTAL_CODING_NT bp"
echo "BWA index action:    $BWA_ACTION"
echo "FAI index action:    $FAI_ACTION"
echo "Reference QC:        $REFERENCE_QC"
echo "Reference CDS QC:    $REFERENCE_CDS_QC"
echo "Index manifest:      $INDEX_MANIFEST"
