#!/bin/bash
# ==========================================================
# author: Meethila Gade
# created on: 31-01-2026-15h-32m
# github: https://github.com/brain-discourse
# copyright: 2026
# description: phase het snp panel
# ==========================================================
#SBATCH -J A2_phase
#SBATCH -o logs/A2_phase.%j.out
#SBATCH -e logs/A2_phase.%j.err
#SBATCH -t 02-00:00:00
#SBATCH -c 2
#SBATCH --mem=80G
module load python/3.12.4
module load samtools
module load gatk

if [[ -n "${EXTRA_PATH:-}" ]]; then export PATH="${EXTRA_PATH}:${PATH}"; fi
which whatshap
whatshap --version

# Project root
PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
DIR="${DIR:-$PROJECT_ROOT}"
MANIFEST="${MANIFEST:-$DIR/config/samples.tsv}"
REF="${REF_FASTA:?Set REF_FASTA to the GRCh38 FASTA used for alignment (with .fai)}"

mkdir -p "$DIR/work/phase" "$DIR/logs"

# One job per case
: "${CASE:?Set CASE to a case label from the manifest, e.g. sbatch --export=ALL,CASE=CASE_A}"
ROW=$(awk -F'\t' -v c="$CASE" 'NR>1 && $2=="bulk" && $7==c{print; exit}' "$MANIFEST")
[[ -n "$ROW" ]] || { echo "[ERROR] No bulk row for CASE=$CASE in $MANIFEST" >&2; exit 2; }

BULK_SAMPLE=$(echo "$ROW" | cut -f1)
BULK_BAM=$(echo "$ROW" | cut -f3)

IN_VCF=$DIR/work/hets/${BULK_SAMPLE}.chr3.hets.tight.vcf.gz
OUT_VCF=$DIR/work/phase/${BULK_SAMPLE}.chr3.hets.tight.phased.vcf.gz
OUT_STATS=$DIR/work/phase/${BULK_SAMPLE}.chr3.hets.tight.phased.stats.tsv

[[ -s "$BULK_BAM" ]] || { echo "[ERROR] Missing BAM: $BULK_BAM" >&2; exit 1; }
[[ -s "$IN_VCF" ]] || { echo "[ERROR] Missing input VCF: $IN_VCF" >&2; exit 1; }
[[ -s "$IN_VCF.tbi" || -s "$IN_VCF.csi" ]] || { echo "[ERROR] Missing VCF index for: $IN_VCF" >&2; exit 1; }
[[ -s "$REF" ]] || { echo "[ERROR] Missing reference fasta: $REF" >&2; exit 1; }
[[ -s "${REF}.fai" ]] || { echo "[ERROR] Missing reference index: ${REF}.fai" >&2; exit 1; }

echo "[DEBUG] ROW=$ROW"
echo "[DEBUG] BULK_SAMPLE=$BULK_SAMPLE"
echo "[DEBUG] BULK_BAM=$BULK_BAM"
echo "[DEBUG] IN_VCF=$IN_VCF"
ls -lh "$BULK_BAM" "$IN_VCF" || exit 1


if [[ ! -s "${BULK_BAM}.bai" ]]; then
  samtools index "$BULK_BAM"
fi

TMP_VCF=$DIR/work/phase/${BULK_SAMPLE}.tmp.phased.vcf

# Phase snps (read-backed; single run, then compress and index)
stdbuf -oL -eL whatshap phase --reference "$REF" -o "$TMP_VCF" "$IN_VCF" "$BULK_BAM"

bgzip -f "$TMP_VCF"
mv -f "${TMP_VCF}.gz" "$OUT_VCF"
tabix -f -p vcf "$OUT_VCF"

whatshap stats --tsv "$OUT_STATS" "$OUT_VCF"

echo "DONE" > $DIR/work/phase/${BULK_SAMPLE}.DONE.ok
